import argparse
parser = argparse.ArgumentParser(description='python msa2vcf.py -i msa.fa ')
parser.add_argument("-i", type=str,help= 'input file')
args = parser.parse_args()
fai=1
faa=''
fab=''
fac=''
fad=''
with open(args.i, 'r+', encoding='UTF-8') as f:
    while 1:
        fas=f.readline()
        if fai==1:
            faa=fas.replace("\n",'').replace(">",'')
        if fai==2:
            fab=fas.replace("\n",'')
        if fai==3:
            fac=fas.replace("\n",'').replace(">",'')
        if fai==4:
            fad=fas.replace("\n",'')
        fai=fai+1
        if not fas:
            break
vcf1 =''
for i,x in enumerate(fab):
    fass = fab[i]+' '+fad[i]+'\n'
    vcf1 = vcf1+ fass
vcf2 =''
bj = vcf1.split('\n')
bjj = []
for line1 in bj:
    line1 = line1.split()
    if len(line1) != 0:
        bjj.append(line1)
for i in bjj:
    if i[0] == i[1] and i[0] != '-' and i[1] != '-' :
        bjj1=i[0]+' '+i[1]+' '+'0'+'\n'
        vcf2 = vcf2+bjj1
    if i[0] != i[1] and i[0] != '-' and i[1] != '-':
        bjj2=i[0]+' '+i[1]+' '+'1'+'\n'
        vcf2 = vcf2 + bjj2
    if i[0] != '-' and i[1] == '-' or i[1] != '-' and i[0] == '-':
        bjj3=i[0]+' '+i[1]+' '+ '3' + '\n'
        vcf2 = vcf2 + bjj3
vcf3 = ''
l2h = vcf2.split('\n')
l2hh = []
for line2 in l2h:
    line2 = line2.split()
    if len(line2) != 0:
        l2hh.append(line2)
l2ha =''
l2hb =''
l2hc =''
for i in l2hh:
    l2t1 =f'{i[0]}'
    l2ha =l2ha+l2t1
    l2t2 =f'{i[1]}'
    l2hb =l2hb+l2t2
    l2t3 =f'{i[2]}'
    l2hc =l2hc+l2t3
vcf3 = l2ha+'\n'+l2hb+'\n'+l2hc+'\n'
vcf4 = vcf3.replace('03','33')
vcf4 = vcf4.replace('13','33')
vcf4 = vcf4.replace('30','34')
vcf4 = vcf4.replace('31','34')
vcf5=''
h2la=''
h2lb=''
h2lc=''
h2l = vcf4.split('\n')
h2la = h2la+ h2l[0]
h2lb = h2lb+ h2l[1]
h2lc = h2lc+ h2l[2]
for i,x in enumerate(h2la):
    h2ll=h2la[i] + " " + h2lb[i]+" "+h2lc[i]+"\n"
    vcf5 = vcf5 + h2ll
wda = 1
wdls =[]
vcf6 = ''
wd = vcf5.split('\n')
for line3 in wd:
    line3 = line3.split()
    if len(line3) != 0:
        wdls.append(line3)
for i in wdls:
    if i[0] != '-' and i[0] != '' :
        wd1 = str(wda)+" "+i[0]+" "+i[1]+" "+i[2]+"\n"
        vcf6 = vcf6+wd1
        wda=wda+1
    if i[0] == '-':
        wd2="BBB"+" "+i[0]+" "+i[1]+" "+i[2]+"\n"
        vcf6 = vcf6 + wd2
vcf7 =''
vcff = []
vcf = vcf6.split('\n')
for line4 in vcf:
    line4 = line4.split()
    if len(line4) != 0:
        vcff.append(line4)
a = ''
b = ''
c = []
for i in vcff:
    if int(i[3]) == int(1):
        tq1 = i[0]+'\t'+i[1]+'\t'+i[2]+'\n'
        vcf7 = vcf7+tq1
    if int(i[3]) == int(3):
        t1 = f'{i[1]}'
        t2 = f'{i[2]}'
        t3 = f'{i[0]}'
        a =a+t1
        b =b+t2
        c.append(t3)
    if int(i[3]) == int(4):
        tq2 = str(c[0]) + '\t' + a + '\t' + b + '\n'
        vcf7 = vcf7 + tq2
        a =''
        b =''
        c =[]
        continue
vcf8 = vcf7.replace('-','')
vcf9 =''
tccw1 = []
tccw = vcf8.split('\n')
for line5 in tccw:
    line5 = line5.split()
    if len(line5) != 0:
        tccw1.append(line5)
for i in tccw1:
    if len(i[1]) !=1 and len(i[2]) !=1:
        sssssssssss = 1
    else:
        tccw2 = f'{faa}'+'\t'+i[0]+'\t'+'.'+'\t'+i[1]+'\t'+i[2]+'\t'+'.'+'\t'+'PASS'+'\t'+'DP=1;AF=1.0'+'\t'+'GT:PL'+'\t'+'1/1:1,22'+'\n'
        vcf9 = vcf9+tccw2

with open(f'{fac}.vcf', 'a', encoding='UTF-8') as fw:
    titou = f'##fileformat=VCFv4.2\n##source=msa2vcf\n##contig=<ID=BJ3>\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{fac}\n'
    fw.write(titou)
    fw.write(vcf9)
print(f'{fac} success')
