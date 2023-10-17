import pysam as ps
import pybedtools as pb
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
group_input = parser.add_argument_group('input (mandatory)')
group_input.add_argument('bamfile_IN', help='BAM file with reads')
group_input.add_argument('bamfile_OUT', help='BAM file with deduplicated reads')
params = parser.parse_args()



print("reading alimnments")
bam_sam = ps.AlignmentFile(params.bamfile_IN,"rb")
bam_sam_out = ps.AlignmentFile(params.bamfile_OUT, "wb", template=bam_sam)
print("reading BEDtool")
bam_bed = pb.BedTool(params.bamfile_IN)
print("converting to BED")
bed = bam_bed.bam_to_bed()
print("creating DataFrame")
bed_pd = bed.to_dataframe()
print("separating ID")
bed_pd['id'],bed_pd['tags'],bed_pd['barcode'],bed_pd['quality'] = bed_pd['name'].str.split("#",3).str
print("splitting strands")
bed_pd_forward = bed_pd[bed_pd['strand'] == "+"]
bed_pd_reverse = bed_pd[bed_pd['strand'] == "-"]
print("grouping")
groups_forward = bed_pd_forward.groupby(['chrom','start','barcode'])
groups_reverse = bed_pd_reverse.groupby(['chrom','end','barcode'])

positive_ids = {}
print("deduping forward strand")
for index,group in groups_forward:
    positive_ids[group.iloc[0]['name']]=True
print("deduping reverse strand")
for index,group in groups_reverse:
    positive_ids[group.iloc[0]['name']] = True
print("### deduplicated entries:", len(positive_ids))
print("saving deduplicated into new BAM")
for map in bam_sam:
    #print(map.query_name)
    #print(map.query_name in positive_ids)
    if(map.query_name in positive_ids):
        bam_sam_out.write(map)

bam_sam_out.close()
bam_sam.close()
