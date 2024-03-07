mkdir deepvariant-survindel2-DEL-30-to-50bp
mkdir deepvariant-survindel2-INS-30-to-50bp
for f in deepvariant-DEL-30-to-50bp/*.vcf.gz ; do bcftools concat -a $f survindel2-DEL-30-to-50bp/`basename $f` -Oz -o deepvariant-survindel2-DEL-30-to-50bp/`basename $f`; done
for f in deepvariant-INS-30-to-50bp/*.vcf.gz ; do bcftools concat -a $f survindel2-INS-30-to-50bp/`basename $f` -Oz -o deepvariant-survindel2-INS-30-to-50bp/`basename $f`; done
