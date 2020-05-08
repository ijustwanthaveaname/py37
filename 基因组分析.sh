# installing MCL (https://www.micans.org/mcl/)
#wget https://www.micans.org/mcl/src/mcl-14-137.tar.gz -P ~/software/
tar zxf ~/software/mcl-14-137.tar.gz
cd mcl-14-137/
./configure --prefix=/opt/biosoft/mcl-14-137/ && make -j 4 && make install
cd .. && rm -rf mcl-14-137/
echo 'PATH=$PATH:/opt/biosoft/mcl-14-137/bin/' >> ~/.bashrc 
source ~/.bashrc 

# installing OrthoMCL (http://orthomcl.org/orthomcl/)
#wget http://orthomcl.org/common/downloads/software/v2.0/orthomclSoftware-v2.0.9.tar.gz -P ~/software/l
tar zxf ~/software/orthomclSoftware-v2.0.9.tar.gz -C /opt/biosoft/
echo 'PATH=$PATH:/opt/biosoft/orthomclSoftware-v2.0.9/bin/' >> ~/.bashrc
source ~/.bashrc

# installing mafft (https://mafft.cbrc.jp/alignment/software/)
#wget https://mafft.cbrc.jp/alignment/software/mafft-7.407-without-extensions-src.tgz -P ~/software
tar zxf ~/software/mafft-7.407-without-extensions-src.tgz
cd mafft-7.407-without-extensions/core/
perl -p -i -e 's#PREFIX =.*#PREFIX = /opt/biosoft/mafft#' Makefile
perl -p -i -e 's#BINDIR =.*#BINDIR = /opt/biosoft/mafft/bin/#' Makefile
make -j 4
make install
cd ../../ && rm -rf mafft-7.407-without-extensions
echo 'PATH=$PATH:/opt/biosoft/mafft/bin/' >> ~/.bashrc
source ~/.bashrc

# installing Gblocks (http://molevol.cmima.csic.es/castresana/Gblocks.html)
#wget http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_Linux64_0.91b.tar.Z -P ~/software/
tar zxf ~/software/Gblocks_0.91b.tar.gz -C /opt/biosoft/
echo 'PATH=$PATH:/opt/biosoft/Gblocks_0.91b/' >> ~/.bashrc
source ~/.bashrc

# installing PAML
#wget http://abacus.gene.ucl.ac.uk/software/paml4.9i.tgz -P ~/software/
tar zxf ~/software/paml4.9i.tgz -C /opt/biosoft/
cd /opt/biosoft/paml4.9i/
rm bin/*
cd src
make -f Makefile
cp baseml basemlg chi2 codeml evolver infinitesites mcmctree pamp yn00 ../bin
echo 'PATH=$PATH:/opt/biosoft/paml4.9i/bin' >> ~/.bashrc
source ~/.bashrc

# installing prottest (https://github.com/ddarriba/prottest3/releases)
#wget https://github.com/ddarriba/prottest3/releases/download/3.4.2-release/prottest-3.4.2-20160508.tar.gz -P ~/software/
tar zxf ~/software/prottest-3.4.2-20160508.tar.gz -C /opt/biosoft/
echo 'export PROTTEST_HOME=/opt/biosoft/prottest-3.4.2' >> ~/.bashrc
source ~/.bashrc

# installing RAxML (https://github.com/stamatak/standard-RAxML | https://cme.h-its.org/exelixis/web/software/raxml/index.html)
#wget https://github.com/stamatak/standard-RAxML/archive/v8.2.12.tar.gz -O ~/software/RAxML-v8.2.12.tar.gz
tar zxf ~/software/RAxML-v8.2.12.tar.gz -C /opt/biosoft/
mv /opt/biosoft/standard-RAxML-8.2.12/ /opt/biosoft/RAxML-8.2.12/
cd /opt/biosoft/RAxML-8.2.12/
make -f Makefile.SSE3.PTHREADS.gcc -j 4
rm *.o
make -f Makefile.AVX.PTHREADS.gcc -j 4
rm *.o
source ~/.bashrc.mpich
make -f Makefile.SSE3.HYBRID.gcc -j 4
rm *.o
make -f Makefile.AVX.HYBRID.gcc -j 4
rm *.o
chmod 755 /opt/biosoft/RAxML-8.2.12/usefulScripts/*
echo 'PATH=$PATH:/opt/biosoft/RAxML-8.2.12/' >> ~/.bashrc
source ~/.bashrc

# installing FigTree (http://tree.bio.ed.ac.uk/software/figtree/ | https://github.com/rambaut/figtree/releases)
#wget https://github.com/rambaut/figtree/releases/download/v1.4.4/FigTree_v1.4.4.tgz -P ~/software
tar zxf ~/software/FigTree_v1.4.4.tgz -C /opt/biosoft/

# installing r8s (https://sourceforge.net/projects/r8s/)
#wget https://sourceforge.net/projects/r8s/files/r8s1.81.tar.gz -P ~/software/
tar zxf ~/software/r8s1.81.tar.gz -C /opt/biosoft/
cd /opt/biosoft/r8s1.81/src
make -j 4
echo 'PATH=$PATH:/opt/biosoft/r8s1.81/src' >> ~/.bashrc
source ~/.bashrc

# installing BEAST2 (https://www.beast2.org/ | https://github.com/CompEvol/beast2)
#wget https://github.com/CompEvol/beast2/releases/download/v2.5.2/BEAST.v2.5.2.Linux.tgz -P ~/software/
tar zxf ~/software/BEAST.v2.5.2.Linux.tgz -C /opt/biosoft/
echo 'PATH=$PATH:/opt/biosoft/beast/bin/' >> ~/.bashrc
source ~/.bashrc

# installing BEAGLE (https://github.com/beagle-dev/beagle-lib)
#wget https://github.com/beagle-dev/beagle-lib/archive/v3.1.2.tar.gz -O ~/software/beagle-lib-3.1.2.tar.gz
tar zxf ~/software/beagle-lib-3.1.2.tar.gz
cd beagle-lib-3.1.2/
./autogen.sh 
./configure --prefix=/opt/biosoft/beagle-lib-3.1.2
make -j 8
make install
echo 'export PKG_CONFIG_PATH=/opt/biosoft/beagle-lib-3.1.2/lib/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=/opt/biosoft/beagle-lib-3.1.2/lib/:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=/opt/biosoft/beagle-lib-3.1.2/include:$C_INCLUDE_PATH' >> ~/.bashrc

# installing Tracer (http://tree.bio.ed.ac.uk/software/tracer/)
#wget https://github.com/beast-dev/tracer/releases/download/v1.7.1/Tracer_v1.7.1.tgz -P ~/software
tar zxf ~/software/Tracer_v1.7.1.tgz -C /opt/biosoft/
chmod 755 /opt/biosoft/Tracer_v1.7.1/bin/tracer 
echo 'PATH=$PATH:/opt/biosoft/Tracer_v1.7.1/bin/' >> ~/.bashrc
source ~/.bashrc


# installing CAFE (https://github.com/hahnlab/CAFE)
#wget https://github.com/hahnlab/CAFE/archive/v4.2.1.tar.gz -O ~/software/CAFE-4.2.1.tar.gz
tar zxf ~/software/CAFE-4.2.1.tar.gz -C /opt/biosoft
cd /opt/biosoft/CAFE-4.2.1
./configure && make -j 4
mkdir bin
cp cafe/caferror.py release/cafe bin/
echo 'PATH=$PATH:/opt/biosoft/CAFE-4.2.1/bin/' >> ~/.bashrc
source ~/.bashrc

# install MCScanX (http://chibba.pgml.uga.edu/mcscan2/)
#wget http://chibba.pgml.uga.edu/mcscan2/MCScanX.zip -P ~/software/
unzip ~/software/MCScanX.zip -d /opt/biosoft/
cd /opt/biosoft/MCScanX/
#make
echo 'PATH=$PATH:/opt/biosoft/MCScanX/' >> ~/.bashrc
source ~/.bashrc

# installing MUMmer (https://github.com/mummer4/mummer)
#wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz -P ~/software
tar zxf ~/software/mummer-4.0.0beta2.tar.gz
cd mummer-4.0.0beta2
./configure --prefix=/opt/biosoft/mummer-4.0.0beta2 && make -j 24 && make install
cd .. && rm -rf mummer-4.0.0beta2
echo 'PATH=$PATH:/opt/biosoft/mummer-4.0.0beta2/bin/' >> ~/.bashrc
source ~/.bashrc

# installing Mauve (http://darlinglab.org/mauve/mauve.html)
#wget http://darlinglab.org/mauve/downloads/mauve_linux_2.4.0.tar.gz -P ~/software/
tar zxf ~/software/mauve_linux_2.4.0.tar.gz -C /opt/biosoft/
echo 'PATH=$PATH:/opt/biosoft/mauve_2.4.0/' >> ~/.bashrc
source ~/.bashrc
####运行脚本
## a. 准备比较基因分析多物种的FASTA序列文件
mkdir -p /home/train/14.genome_comparison/a.preparing_data
cd /home/train/14.genome_comparison/a.preparing_data

# 1. 通过NCBI的Genome数据库(https://www.ncbi.nlm.nih.gov/genome)查找目标物种的信息
# 手工生成文件source.txt，该文件前4列为：物种缩写名（推荐为属名前2个字母+种名前3个字母）、双名法拉丁文物种名、NCBI的Genome数据库网址、从NCBI下载的数据文件前缀。
# 该文件前第一列用于后续的数据文件名前缀、第二列用于后续将结果中的物种名缩写转换成物种名全称、第三列用于数据的批量化下载。此外，手动编辑文件source.txt时，可以生成更多的列，记录其它基因组信息。
echo -e "parub\tPaxillus rubicundulus\thttps://www.ncbi.nlm.nih.gov/genome/35925
sccit\tScleroderma citrinum\thttps://www.ncbi.nlm.nih.gov/genome/18213
laame\tLaccaria amethystina\thttps://www.ncbi.nlm.nih.gov/genome/17383
plost\tPleurotus ostreatus\thttps://www.ncbi.nlm.nih.gov/genome/909?genome_assembly_id=202364
lasul\tLaetiporus sulphureus\thttps://www.ncbi.nlm.nih.gov/genome/17416
phgig\tPhlebiopsis gigantea\thttps://www.ncbi.nlm.nih.gov/genome/18214" > source.txt

# 2. 批量化下载数据：基因组序列和基因结构注释文件
#perl -e 'while (<>) { chomp; @_ = split /\t/; $curl = `curl $_[2]`; if ($curl =~ m/Download genome annotation in.*href=\"(.*)\">GFF/) { $file = $1; $file =~ s/_genomic.gff.gz//; $prefix = $file; $prefix =~ s/.*\///; print "wget $file\_genomic.gff.gz -O $_[0].genomic.gff.gz\nwget $file\_genomic.fna.gz -O $_[0].genomic.fna.gz\n"; } }' source.txt > command.download.list 2> /dev/null
#ParaFly -c command.download.list -CPU 4
tar zxf /home/train/00.incipient_data/data_for_genome_comparison/genomics_data_from_NCBI.tar.gz

# 3. 对NCBI的GFF文件进行格式修正。默认下NCBI的GFF格式不太正确，且其基因ID编号比较乱。
# 对NCBI的GFF3文件进行整理，仅保留含有mRNA信息的基因，修正尾部的CDS（有些NCBI的注释结果中最后一个CDS是不包含终止密码子的），修正同一条链上存在重叠的基因情况，并对基因id进行重命名。
# 最后，得到所有基因的Protein和CDS序列。若一个基因存在多个可变剪接，则仅保留CDS最长的转录本信息。
for i in `cut -f 1 source.txt`
do
    echo "gzip -dc $i.genomic.gff.gz > $i.genome.gff; gzip -dc $i.genomic.fna.gz > $i.genome.fasta; gff3_clear.pl --prefix $i $i.genome.gff > $i.genome.gff3;  GFF3Clear --gene_prefix $i --genome $i.genome.fasta $i.genome.gff3 > $i.geneModels.gff3; gff3ToGtf.pl $i.genome.fasta $i.geneModels.gff3 > $i.geneModels.gtf; eukaryotic_gene_model_statistics.pl $i.geneModels.gtf $i.genome.fasta $i > $i.statistics"
done > command.geneModels.list
ParaFly -c command.geneModels.list -CPU 6


## b. 使用OrthoMCL进行同源基因聚类分析
mkdir -p /home/train/14.genome_comparison/b.orthomcl
cd /home/train/14.genome_comparison/b.orthomcl

# 1. 在当前工作目录中创建配置文件。配置数据库信息
perl -p -e 's/:3307//; s/^dbLogin=.*/dbLogin=train/; s/^dbPassword=.*/dbPassword=123456/' /opt/biosoft/orthomclSoftware-v2.0.9/doc/OrthoMCLEngine/Main/orthomcl.config.template > orthomcl.config

# 2. 安装 OrthoMCL 运行所需要的 mysql 数据库的表
echo "DROP DATABASE orthomcl" | mysql -utrain -p123456
echo "CREATE DATABASE IF NOT EXISTS orthomcl" | mysql -utrain -p123456
orthomclInstallSchema orthomcl.config 

# 3. 准备 OrthoMCL 的输入文件
# 输入文件是包含多个 Fasta 文件的文件夹，每个 Fasta 文件是一个物种的蛋白质组序列。每个 Fasta 文件的序列名必须满足格式, 示例: >ncra|NCU02040 。用 '|' 分开，前者代表物种名，推荐和 Fasta 文件名前缀一致，后者代表独一无二蛋白质ID。 可以使用 orthomclAdjustFasta 对一般的 Fasta 文件进行处理，使之和 OrthoMCL 软件兼容。
mkdir compliantFasta
cd compliantFasta
for i in `cut -f 1 ../../a.preparing_data/source.txt`
do
    echo "orthomclAdjustFasta $i ../../a.preparing_data/$i.pep.fasta 1"
done | sh
cd ..
perl -p -i -e 's/\*$//; s/\*/X/g' compliantFasta/*.fasta
mkdir compliantFasta_CDS
cd compliantFasta_CDS
for i in `cut -f 1 ../../a.preparing_data/source.txt`
do
    echo "orthomclAdjustFasta $i ../../a.preparing_data/$i.CDS.fasta 1"
done | sh
cd ..

# 4. 过滤低质量序列。允许的最短的 protien 长度是 30，终止密码子最大比例为 20% 。
orthomclFilterFasta compliantFasta/ 30 20
# 该命令只能过滤低质量序列。而输入文件中最好还需要过滤掉可变剪切。

# 5. 对过滤后的结果 goodProteins.fasta 进行 all-vs-all BLAST
diamond makedb --in goodProteins.fasta -d goodProteins
diamond blastp --db goodProteins --query goodProteins.fasta --out diamond.xml --outfmt 5 --sensitive --max-target-seqs 500 --evalue 1e-5 --id 20 --tmpdir /dev/shm
parsing_blast_result.pl --no-header --max-hit-num 500 --evalue 1e-9 --identity 0.1 --CIP 0.3 --subject-coverage 0.5 --query-coverage 0.5 diamond.xml > diamond.tab
#gzip -dc /home/train/00.incipient_data/data_for_genome_comparison/diamond.tab.gz > diamond.tab

# 6. 对 Blast 的结果进行处理，得到序列两两之间的相似性信息，以利于导入到 mysql 
orthomclBlastParser diamond.tab compliantFasta/ > similarSequences.txt
# similarSequences.txt 文件内容有 8 列：query_id, subject_id, query_taxon, subject_taxon, evalue_mant, evalue_exp, percent_ident 和 percent_match 。

# 7. 将 similarSequences.txt 导入到 mysql 数据库
orthomclLoadBlast orthomcl.config similarSequences.txt

# 8. 寻找相似序列对：在物种间中是相互最佳匹配的对（orthologs对）；在物种内是相互最佳匹配并优于种间匹配的对（in-paralogs对）；前两者结合得到的co-orthologs对。
orthomclPairs orthomcl.config orthomcl_pairs.log cleanup=no

# 9. 将找到的相似序列对从 mysql 中导出
orthomclDumpPairsFiles orthomcl.config

# 10. 使用 mcl 对 pairs 进行聚类（Ortholog Cluster Groups）并对类进行编号
mcl mclInput --abc -I 1.5 -o mclOutput
orthomclMclToGroups OCG 1 < mclOutput > groups.txt
perl -e 'open IN, $ARGV[0]; @num = <IN>; $num = @num; $len = length($num); foreach (@num) { m/OCG(\d+)/; $id = 0 x ($len - length($1)); s/OCG/OCG$id/; print;} ' groups.txt > aa; mv aa groups.txt

# 11. 对 groups.txt 进行同源基因数量的统计
orthomcl_genes_number_stats.pl groups.txt compliantFasta > genes_number_stats.txt

# 12. 进行旁系同源基因分析
orthomcl_get_outParalog.pl groups.txt similarSequences.txt
cd ..


## c. 使用单拷贝同源基因构建物种树
mkdir -p /home/train/14.genome_comparison/c.species_tree
cd /home/train/14.genome_comparison/c.species_tree

# 1. 根据 orthomcl 结果提取单拷贝同源基因
orthomcl_extract_ortholog_seqs.pl --out_directory orthologGroups_CDS --species_ratio 1 --single_copy_species_ratio 0.5 --copy_num 10 --max_seq_length 6000 ../b.orthomcl/groups.txt ../b.orthomcl/compliantFasta_CDS
orthomcl_extract_ortholog_seqs.pl --out_directory orthologGroups_Protein --species_ratio 1 --single_copy_species_ratio 0.5 --copy_num 10 --max_seq_length 6000 --compliantFasta_dir_for_calculating_seq_length  ../b.orthomcl/compliantFasta_CDS ../b.orthomcl/groups.txt ../b.orthomcl/compliantFasta
perl -p -i -e 's/\*$//; s/\*/X/g' orthologGroups_Protein/*.fasta
ls orthologGroups_Protein/ | perl -pe 's/.fasta//' > orthologGroups.txt
perl -e 'while (<>) { $num ++; print if ($num % 10) == 1; }' orthologGroups.txt > aa; mv aa orthologGroups.txt

# 2. 对单拷贝同源进行进行多序列比对
for i in `cat orthologGroups.txt`
do
    echo "linsi orthologGroups_Protein/$i.fasta > orthologGroups_Protein/$i.fasta.align"
done > command.mafft.list
ParaFly -c command.mafft.list -CPU 8

# 3. 将Protein序列比对结果转换为Codon序列比对结果
for i in `cat orthologGroups.txt`
do
    echo "proteinAlignment2CDSAlignemnt.pl orthologGroups_Protein/$i.fasta.align orthologGroups_CDS/$i.fasta > orthologGroups_CDS/$i.fasta.align"
done > command.proteinAlignment2CDSAlignemnt.list
ParaFly -c command.proteinAlignment2CDSAlignemnt.list -CPU 8

# 4. 对Codon序列比对结果进行保守区块提取
for i in `cat orthologGroups.txt`
do
    echo "Gblocks orthologGroups_CDS/$i.fasta.align -t=c; if [ -f orthologGroups_CDS/$i.fasta.align-gb ]; then echo $i completed; fi"
    echo "Gblocks orthologGroups_Protein/$i.fasta.align -t=p; if [ -f orthologGroups_Protein/$i.fasta.align-gb ]; then echo $i completed; fi"
done > command.Gblocks.list
ParaFly -c command.Gblocks.list -CPU 8

# 5. 将各个同源基因家族的多序列比对结果转换为Phylip格式
for i in `cat orthologGroups.txt`
do
    perl -p -e 's/(^>\w+).*/$1/; s/ //g' orthologGroups_CDS/$i.fasta.align-gb > orthologGroups_CDS/$i.fasta.align-gb.fasta
done
for i in `cat orthologGroups.txt`
do
    echo "/opt/biosoft/RAxML-8.2.12/usefulScripts/convertFasta2Phylip.sh orthologGroups_CDS/$i.fasta.align-gb.fasta > orthologGroups_CDS/$i.phy"
done > command.convertFasta2Phylip.list
ParaFly -c command.convertFasta2Phylip.list -CPU 8

# 6. 整合所有单拷贝同源基因的Codon序列比对结果，并使用RAxML构建物系统发育树
mkdir RAxML
perl -e 'while (<orthologGroups_CDS/*.align-gb>) { open IN, $_ or die $!; while (<IN>) { if (m/^>(\w+)/) { $seq_id = $1; } else { s/\s+//g; $seq{$seq_id} .= $_; } } } foreach (sort keys %seq) { print ">$_\n$seq{$_}\n"; }' > RAxML/allSingleCopyOrthologsAlign.Codon.fasta
cd RAxML
/opt/biosoft/RAxML-8.2.12/usefulScripts/convertFasta2Phylip.sh allSingleCopyOrthologsAlign.Codon.fasta > allSingleCopyOrthologsAlign.Codon.phy
/opt/biosoft/RAxML-8.2.12/raxmlHPC-PTHREADS-SSE3 -f a -x 12345 -p 12345 -# 100 -m GTRGAMMA -s allSingleCopyOrthologsAlign.Codon.phy -n out_codon -T 8
# FigTree画树，并输出树文件信息
java -jar /opt/biosoft/FigTree_v1.4.4/lib/figtree.jar RAxML_bipartitions.out_codon
cd ../

# 7. 可能有部分单拷贝同源基因不利于构建系统发育树，将这些基因挑出来去除掉后，有利于构建更准确的物种树
# 使用PAML软件的baseml对各个单拷贝同源基因进行枝长分析，再与标准树进行比较，剔除异常的基因。
mkdir paml_baseml
cd paml_baseml
# 7.1 根据RAxML的结果提取系统发育树的拓扑结构
echo "6 1" > input.trees
echo "(((parub,sccit),(laame,plost)),(phgig,lasul));" >> input.trees
# 7.2 分别对各个单拷贝同源基因的密码子比对结果用codeml进行系统发育分析，得到每个各个同源基因的树信息
mkdir OCG_trees
for i in `cat ../orthologGroups.txt`
do
    echo "calculating_branchLength_by_baseml.pl ../orthologGroups_CDS/$i.phy input.trees > OCG_trees/$i.tree"
done > command.calculating_branchLength_by_baseml.list
ParaFly -c command.calculating_branchLength_by_baseml.list -CPU 8
# 7.3 分析每个各个同源基因树枝长和RAxML树的平均偏差百分比：单拷贝同源基因某分枝的偏差百分比 = (绝对值{单拷贝同源基因树某分枝的枝长 - RAxML树中该分枝的枝长} / RAxML树中该分枝的枝长) * 100%；然后计算所有分枝偏差百分比的平均值。
echo "(((parub:0.2669744817942059,sccit:0.33423852353852823):0.2746796113280112,(laame:0.4546314656446836,plost:0.48573483415848673):0.11098577281827804):0.09014410998935125,(phgig:0.46315428814836046,lasul:0.3948929446629426):0.09014410998935123);" > standard.tree
for i in `cat ../orthologGroups.txt`
do
    echo "calculating_branchLength_bias_percentage_of_two_trees.pl standard.tree OCG_trees/$i.tree > OCG_trees/$i.bias"
done > command.calculating_branchLength_bias_percentage_of_two_trees.list
ParaFly -c command.calculating_branchLength_bias_percentage_of_two_trees.list -CPU 8
perl -e 'while (<OCG_trees/*.bias>) { my $ocg = $1 if m/(\w+).bias/; open IN, $_; $_ = <IN>; s/\%//g; print "$ocg\t$_" if $_; }' > aa;
remove_extremum_through_normal_distribution.pl --data_column 4 --confidence_interval 2 aa > bb
remove_extremum_through_normal_distribution.pl --data_column 4 --confidence_interval 2 bb > cc
cut -f 1 cc > keep_OCG.txt
cd ..

# 8. 再次使用RAxML构建系统发育树
mkdir RAxML_again
perl -e 'open IN, "paml_baseml/keep_OCG.txt"; while (<IN>) { chomp; push @ocg, $_; } close IN; foreach (@ocg) { open IN, "orthologGroups_CDS/$_.fasta.align-gb"; while (<IN>) { if (m/^>(\w+)/) { $seq_id = $1; } else { s/\s+//g; $seq{$seq_id} .= $_; } } } foreach (sort keys %seq) { print ">$_\n$seq{$_}\n"; }' > RAxML_again/allSingleCopyOrthologsAlign.Codon.fasta
perl -e 'open IN, "paml_baseml/keep_OCG.txt"; while (<IN>) { chomp; push @ocg, $_; } close IN; foreach (@ocg) { open IN, "orthologGroups_Protein/$_.fasta.align-gb"; while (<IN>) { if (m/^>(\w+)/) { $seq_id = $1; } else { s/\s+//g; $seq{$seq_id} .= $_; } } } foreach (sort keys %seq) { print ">$_\n$seq{$_}\n"; }' > RAxML_again/allSingleCopyOrthologsAlign.Protein.fasta
cd RAxML_again
/opt/biosoft/RAxML-8.2.12/usefulScripts/convertFasta2Phylip.sh allSingleCopyOrthologsAlign.Codon.fasta > allSingleCopyOrthologsAlign.Codon.phy
/opt/biosoft/RAxML-8.2.12/usefulScripts/convertFasta2Phylip.sh allSingleCopyOrthologsAlign.Protein.fasta > allSingleCopyOrthologsAlign.Protein.phy
condon_alignment_partition.pl --out-prefix allSingleCopyOrthologsAlign_Codon allSingleCopyOrthologsAlign.Codon.phy
raxmlHPC-PTHREADS-SSE3 -f a -x 12345 -p 12345 -# 1000 -m GTRGAMMA -s allSingleCopyOrthologsAlign.Codon.phy -n out_codon -T 8
#/opt/sysoft/mpich2-1.5/bin/mpirun -np 20 /opt/biosoft/RAxML-8.2.12/raxmlHPC-HYBRID-SSE3 -f a -x 12345 -p 12345 -# 1000 -m PROTGAMMALGX -s allSingleCopyOrthologsAlign.Protein.phy -n out_protein -T 8
#/opt/sysoft/mpich2-1.5/bin/mpirun -np 20 /opt/biosoft/RAxML-8.2.12/raxmlHPC-HYBRID-SSE3 -f a -x 12345 -p 12345 -# 1000 -m GTRGAMMAX -s allSingleCopyOrthologsAlign.Codon.phy -n out_codon -T 8
#/opt/sysoft/mpich2-1.5/bin/mpirun -np 20 /opt/biosoft/RAxML-8.2.12/raxmlHPC-HYBRID-SSE3 -f a -x 12345 -p 12345 -# 1000 -m GTRGAMMAX -s allSingleCopyOrthologsAlign_Codon.4.phy -n out_FF4 -T 8
#/opt/sysoft/mpich2-1.5/bin/mpirun -np 20 /opt/biosoft/RAxML-8.2.12/raxmlHPC-HYBRID-SSE3 -f a -x 12345 -p 12345 -# 1000 -m GTRGAMMAX -s allSingleCopyOrthologsAlign_Codon.3.phy -n out_codon3 -T 8
#/opt/sysoft/mpich2-1.5/bin/mpirun -np 20 /opt/biosoft/RAxML-8.2.12/raxmlHPC-HYBRID-SSE3 -f a -x 12345 -p 12345 -# 1000 -m GTRGAMMAX -s allSingleCopyOrthologsAlign_Codon.12.phy -n out_codon12 -T 8
# 将树文件中的物种名缩写换成全称
cp RAxML_bipartitions.out_codon tree_abbr.RAxML
cp tree_abbr.RAxML tree_fullName.RAxML
cut -f 1,2 ../../a.preparing_data/source.txt | perl -p -e 's/\t/\t\"/; s/$/\"/;' | perl -p -e 's#\s+#/#; s#^#perl -p -i -e "s/#; s#\n$#/" tree_fullName.RAxML\n#;' | perl -p -e 's#/\"#/\\\"#; s#\"/#\\\"/#;' | sh
# FigTree画树
java -jar /opt/biosoft/FigTree_v1.4.4/lib/figtree.jar tree_fullName.RAxML
cd ../../


# d_1. 使用PAML软件的mcmctree进行物种分歧时间计算
mkdir -p /home/train/14.genome_comparison/d.divergence_time_mcmctree
cd /home/train/14.genome_comparison/d.divergence_time_mcmctree
# (1) 准备Newick格式的树文件，根据RAxML的结果提取系统发育树的拓扑结构
echo "6 1" > input.trees
echo "(((laame,plost),(parub,sccit)'>0.76<1.08')'>2.29<2.66',(lasul,phgig));" >> input.trees
# Paxillus rubicundulus <=> Scleroderma citrinum        89MYA 76-108MYA
# Paxillus rubicundulus <=> Laccaria amethystina	246MYA 229-266MYA
# (2) 准备多序列比对Phylip格式输入文件，可以综合三个位点核酸序列的多序列比对文件
cat ../c.species_tree/RAxML_again/allSingleCopyOrthologsAlign_Codon.[123].phy | perl -pe 's/\s+/  /' > input.txt
# (3) 准备paml mcmctree配置文件
perl -p -e 's/mtCDNApri123/input/; s/mtCDNApri/input/; s/<1.0/<8.0/; s/^(\s*ndata =) .*/$1 3/; s/usedata = .*/usedata = 3/; s/model = .*/model = 7/; s/alpha = .*/alpha = 0.5/; s/ncatG = .*/ncatG = 5/; s/cleandata = .*/cleandata = 0/; s/burnin = .*/burnin = 400000/; s/sampfreq = .*/sampfreq = 100/; s/nsample = .*/nsample = 10000/;' /opt/biosoft/paml4.9i/examples/DatingSoftBound/mcmctree.ctl > mcmctree.ctl
# 推荐加大取样量和burnin数量。为了让程序运行更快，上面的命令减少了10倍。
#perl -p -e 's/mtCDNApri123/input/; s/mtCDNApri/input/; s/<1.0/<8.0/; s/^(\s*ndata =) .*/$1 3/; s/usedata = .*/usedata = 3/; s/model = .*/model = 7/; s/alpha = .*/alpha = 0.5/; s/ncatG = .*/ncatG = 5/; s/cleandata = .*/cleandata = 0/; s/burnin = .*/burnin = 4000000/; s/sampfreq = .*/sampfreq = 100/; s/nsample = .*/nsample = 100000/;' /opt/biosoft/paml4.9i/examples/DatingSoftBound/mcmctree.ctl > mcmctree.ctl
# (4) 运行PAML软件的mcmctree (若数据量较大，则每当程序调用baseml时，按ctrl + c终止，再手动并行化运行baseml)
mcmctree mcmctree.ctl
cp out.BV in.BV
# (5) 手动并行化运行baseml
#echo "mkdir tmp0001; cd tmp0001; ln -s ../input* ../tmp0001* ./; baseml tmp0001.ctl;
#mkdir tmp0002; cd tmp0002; ln -s ../input* ../tmp0002* ./; baseml tmp0002.ctl;
#mkdir tmp0003; cd tmp0003; ln -s ../input* ../tmp0003* ./; baseml tmp0003.ctl;" > command.baseml.list
#ParaFly -c command.baseml.list -CPU 3
#cat tmp0001/rst2 tmp0002/rst2 tmp0003/rst2 > in.BV
# (6) 使用mcmctree进行approximate likelihood分析
perl -p -i -e 's/usedata = .*/usedata = 2    \* 0:/;' mcmctree.ctl
mkdir run01 run02
cp input.txt input.trees mcmctree.ctl in.BV run01
cp input.txt input.trees mcmctree.ctl in.BV run02
echo 'cd run01; mcmctree mcmctree.ctl &> mcmctree.log
cd run02; mcmctree mcmctree.ctl &> mcmctree.log' > command.mcmctree.list
ParaFly -c command.mcmctree.list -CPU 2
# (7) 比较两次运行的MCMC树，若差异较小，则认可其结果
perl -n -e 'my $out; while (s/(.*?(\d\.\d+))//) { my $info = $1; my $value = $2; my $new = $value * 100; $info =~ s/$value/$new/;  $out .= $info; } $out .= $_; print $out' run01/FigTree.tre > tree01.nex
perl -n -e 'my $out; while (s/(.*?(\d\.\d+))//) { my $info = $1; my $value = $2; my $new = $value * 100; $info =~ s/$value/$new/;  $out .= $info; } $out .= $_; print $out' run02/FigTree.tre > tree02.nex
perl -e 'while (<>) { if (s/\s*UTREE.*?=\s*//) { s/\s*\[.*?\]//g; print; } }' tree01.nex > tree01.txt
perl -e 'while (<>) { if (s/\s*UTREE.*?=\s*//) { s/\s*\[.*?\]//g; print; } }' tree02.nex > tree02.txt
calculating_branchLength_bias_percentage_of_two_trees.pl --no_normalization_of_total_branch_length tree01.txt tree02.txt > bias_of_2runs.txt
# 0.41%	0.09%	0.25%
rm tree*
perl -n -e 'my $out; while (s/(.*?(\d\.\d+))//) { my $info = $1; my $value = $2; my $new = $value * 100; $info =~ s/$value/$new/;  $out .= $info; } $out .= $_; print $out' run01/FigTree.tre > tree_abbr.mcmctree
# (8) 将树信息中的简称换成属名和种名全称
cp tree_abbr.mcmctree tree_fullName.mcmctree
cut -f 1,2 ../a.preparing_data/source.txt | perl -p -e 's/\t/\t\"/; s/$/\"/;' | perl -p -e 's#\s+#/#; s#^#perl -p -i -e "s/#; s#\n$#/" tree_fullName.mcmctree\n#;' | perl -p -e 's#/\"#/\\\"#; s#\"/#\\\"/#;' | sh
java -jar /opt/biosoft/FigTree_v1.4.4/lib/figtree.jar tree_fullName.mcmctree
cd ..

## d_2. 使用r8s进行分子钟分析
mkdir -p /home/train/14.genome_comparison/d.divergence_time_r8s
cd /home/train/14.genome_comparison/d.divergence_time_r8s
# 准备r8s的配置文件
echo '#NEXUS
begin trees;
tree tree_1 = [&R] (((laame:0.4619521503678744,plost:0.49080938106855587):0.11126987564556445,(parub:0.2701860446010115,sccit:0.3406314112205876):0.27761850096484975):0.08760670524429123,(lasul:0.3993834859088801,phgig:0.4654508813898232):0.08760670524429123);
end;
begin r8s;
blformat lengths=persite nsites=429963 ulrametric=no;
MRCA BOL sccit parub;
MRCA AGC laame plost parub sccit;
fixage taxon=BOL age=89;
fixage taxon=AGC age=246;' > r8s_in.txt
# 设置多个smoothing参数值进行分歧时间计算，和校准点时刻进行比较，计算其差异大小和差异率。
for ((i=-9;i<=9;i=i+1))
do
    echo "mkdir tmp_$i; cp r8s_in.txt ./tmp_$i/; echo \"divtime method=PL algorithm=TN crossv=yes fossilfixed=yes cvstart=$i cvinc=0.1 cvnum=10;\" >> tmp_$i/r8s_in.txt; r8s -b -f tmp_$i/r8s_in.txt > tmp_$i/r8s_out.txt; echo \"tmp_$i OK\";"
done > command.r8s_crossValidation.list
ParaFly -c command.r8s_crossValidation.list -CPU 19
# 选取差异大小和差异率最小的smoothing
for ((i=-9;i<=5;i=i+1))
do
    grep -P "^\s+-?\d+" tmp_$i/r8s_out.txt
done | sort -k3 -k4 -n | grep Good | head -n 1 | perl -pe 's/\s+/\t/g' | cut -f 3 > best_smoothing_value.txt
# Best smooth value is 
# smooth    Fract_Error    Raw_Error
# 1e-08	    0.3316         51.7657
# Fract_Error表示估算的分歧时间与预设时间之间的差异比例。该值越小越好，若较大，则说明预设时间不准确。
rm command* tmp* -rf
# 设置smoothing参数值
BESTSMOOTH=`cat best_smoothing_value.txt`
echo "set smoothing=$BESTSMOOTH;
divtime method=PL algorithm=TN;
showage;
describe plot=cladogram;
describe plot=chrono_description;
end;" >> r8s_in.txt
# 进行分歧时间计算
r8s -b -f r8s_in.txt > r8s_out.txt
# 将树信息中的简称换成属名和种名全称
perl -e 'while (<>) { if (m/tree tree_1 = (.*);/) { my $tree = $1; $tree =~ s/(AGC|BOL)//g; print "$tree\n"; } }' r8s_out.txt > tree_abbr.r8s
cp tree_abbr.r8s tree_fullName.r8s
cut -f 1,2 ../a.preparing_data/source.txt | perl -p -e 's/\t/\t\"/; s/$/\"/;' | perl -p -e 's#\s+#/#; s#^#perl -p -i -e "s/#; s#\n$#/" tree_fullName.r8s\n#;' | perl -p -e 's#/\"#/\\\"#; s#\"/#\\\"/#;' | sh
java -jar /opt/biosoft/FigTree_v1.4.4/lib/figtree.jar tree_fullName.r8s
cd ..


## d_3. 使用BEAST进行份子钟分析
mkdir -p /home/train/14.genome_comparison/d.divergence_time_BEAST2
cd /home/train/14.genome_comparison/d.divergence_time_BEAST2
phy2nex.pl --data-type DNA ../c.species_tree/RAxML_again/allSingleCopyOrthologsAlign.Codon.phy > input.nex
# 使用beauti命令制作BEAST输入文件input.xml
# Site Model: Substitition Rate estimate 勾上; Gamma Category Count 设置为 4； Subst Model设置为GTR，Rate CT estimate勾上
# Clock Model: Relaxed Clock Log Normal
# Priors: Tree.t => Calibrated Yule Model; birthRateY.t => Gamma(Alpha 0.001, Beta 1000); gammaShape.s => Gamma(Alpha 0.001, Beta 1000)
# BOL(monophyletic): sccit parub	89MYA 76-108MYA	LN(4.485, 0.079)
# AGC(monophyletic): laame plost sccit parub	246MYA 229-266MYA	LN(5.505, 0.037)
# OUT(monophyletic): lasul phgig
# MCMC: Chain Length => 10000000; Store Every => 5000; 
# 最后另存为input.xml
cp /home/train/00.incipient_data/data_for_genome_comparison/input.xml ./
beast -beagle -beagle_CPU -beagle_SSE -beagle_double -threads 8 -instances 8 input.xml &> beast2.log
tracer input.log
treeannotator -burnin 20 input.trees tree_abbr.BEAST2
cp tree_abbr.BEAST2 tree_fullName.BEAST2
cut -f 1,2 ../a.preparing_data/source.txt | perl -p -e 's/\t/\t\"/; s/$/\"/;' | perl -p -e 's#\s+#/#; s#^#perl -p -i -e "s/#; s#\n$#/" tree_fullName.BEAST2\n#;' | perl -p -e 's#/\"#/\\\"#; s#\"/#\\\"/#;' | sh
java -jar /opt/biosoft/FigTree_v1.4.4/lib/figtree.jar tree_fullName.BEAST2
cd ..


## e. 使用CAFE进行基因家族扩张分析
mkdir -p /home/train/14.genome_comparison/e.cafe
cd /home/train/14.genome_comparison/e.cafe

# 根据 orthomcl 结果得到基因家族的数量信息
orthomcl_extract_ortholog_seqs.pl --out_directory orthologGroups_CDS --out_tab_for_CAFE orthomcl2cafe.tab --species_ratio 0.4 --single_copy_species_ratio 0.0 --copy_num 1000 --max_seq_length 100000 ../b.orthomcl/groups.txt ../b.orthomcl/compliantFasta_CDS
rm orthologGroups_CDS -rf
# 提取至少在3个species中存在该 gene family 的OCGs 。

echo '#!/opt/biosoft/CAFE-4.2.1/bin/cafe
version
date

load -i orthomcl2cafe.tab -t 8 -p 0.01
tree (((laame:191.8242,plost:191.8242):44.9287,(parub:107.055,sccit:107.055):129.6979):13.6525,(lasul:187.286,phgig:187.286):63.1194);

lambda -s
report out' > cafe_command

chmod 755 cafe_command
caferror.py -i cafe_command
parsing_cafeOut.pl caferror_1/cafe_final_report.cafe
cd ..


## f. 使用PAML软件的codeml进行正选择基因分析
mkdir -p /home/train/14.genome_comparison/f.PSG_analysis
cd /home/train/14.genome_comparison/f.PSG_analysis
# (1) 准备同源基因的密码子比对结果Phylip格式文件及其sub tree文件
# 选取至少在4个物种中出现的同源基因
orthomcl_extract_ortholog_seqs.pl --out_directory orthologGroups_CDS --species_ratio 0.6 --single_copy_species_ratio 0.0 --copy_num 1000 --max_seq_length 100000 ../b.orthomcl/groups.txt ../b.orthomcl/compliantFasta_CDS
orthomcl_extract_ortholog_seqs.pl --out_directory orthologGroups_Protein --species_ratio 0.6 --single_copy_species_ratio 0.0 --copy_num 1000 --max_seq_length 100000 --compliantFasta_dir_for_calculating_seq_length ../b.orthomcl/compliantFasta_CDS/ ../b.orthomcl/groups.txt ../b.orthomcl/compliantFasta
perl -p -i -e 's/\*/X/g' orthologGroups_Protein/*.fasta
ls orthologGroups_Protein/ | perl -pe 's/.fasta//' > orthologGroups.txt
# 对同源基因的蛋白序列进行多序列比对
for i in `cat orthologGroups.txt`
do
    echo "linsi orthologGroups_Protein/$i.fasta > orthologGroups_Protein/$i.fasta.align"
done > command.mafft.list
ParaFly -c command.mafft.list -CPU 8
# tar zxf /home/train/00.incipient_data/data_for_genome_comparison/PSG_orthologGroups_Protein.tar.gz
# 将蛋白序列比对结果转换为密码子序列比对结果
for i in `cat orthologGroups.txt`
do
    echo "proteinAlignment2CDSAlignemnt.pl orthologGroups_Protein/$i.fasta.align orthologGroups_CDS/$i.fasta > orthologGroups_CDS/$i.fasta.align"
done > command.proteinAlignment2CDSAlignemnt.list
ParaFly -c command.proteinAlignment2CDSAlignemnt.list -CPU 8
perl -p -i -e 's/^(>\w+).*/$1/' orthologGroups_CDS/*.fasta.align
# 分别对每个同源基因进行分析，得到其Phylip格式的多序列比对结果和sub tree文件
echo "(((laame,plost),(parub,sccit)),(lasul,phgig));" >> tree.txt
for i in `cat orthologGroups.txt`
do
    i=${i/*\//}
    i=${i/.fasta.align/}
    echo "/opt/biosoft/RAxML-8.2.12/usefulScripts/convertFasta2Phylip.sh orthologGroups_CDS/$i.fasta.align > orthologGroups_CDS/$i.phy; extract_sub_tree.pl --no_branch_length --no_root orthologGroups_CDS/$i.fasta.align tree.txt > orthologGroups_CDS/$i.tree"
done > command.Phylip_and_subTree.list
ParaFly -c command.Phylip_and_subTree.list -CPU 8
# (2) 通过对同源基因序列之间两两计算dn/ds方法，筛选候选的正选择基因
# 使用YN00算法计算dn/ds
mkdir dnds_yn00
for i in `cat orthologGroups.txt`
do
    echo "calculating_omega_by_yn00.pl --omega_for_PSG 1.0 orthologGroups_CDS/$i.fasta.align > dnds_yn00/$i.txt 2> dnds_yn00/$i.stats"
done > command.calculating_omega_by_yn00.list
ParaFly -c command.calculating_omega_by_yn00.list -CPU 8
ls orthologGroups_CDS/*.fasta.align.PSGyes | perl -ne 'print "$1\n" if m/(OCG\d+)/' > PSG_yn00.list
# 使用ML方法计算dn/ds
# tar zxf /home/train/00.incipient_data/data_for_genome_comparison/PSG_orthologGroups_CDS.tar.gz
mkdir dnds_ML
for i in `cat orthologGroups.txt`
do
    echo "calculating_omega_by_codeml.pl --omega_for_PSG 1.0 orthologGroups_CDS/$i.phy orthologGroups_CDS/$i.tree > dnds_ML/$i.txt 2> dnds_ML/$i.stats"
done > command.calculating_omega_by_ML.list
ParaFly -c command.calculating_omega_by_ML.list -CPU 8
ls orthologGroups_CDS/*.phy.PSGyes | perl -ne 'print "$1\n" if m/(OCG\d+)/' > PSG_ML.list
# 确定待分析的目标分枝
perl -e 'while (<>) { @_ = m/(\w+)/g; } foreach (@_) { print "$_:$_\n"; }' tree.txt > branches_outerNode.txt
echo "Agaricales:laame,plost
Boletales:parub,sccit
Polyporales:lasul,phgig" > branches_innerNode.txt
# 分析各分枝的候选正选择基因
for x in `cat branches_outerNode.txt branches_innerNode.txt`
do
    y=$x
    x=${x/:*/}
    y=${y/*:/}

    mkdir BS_$x; cd BS_$x; mkdir dnds_yn00
    for i in `cat ../orthologGroups.txt`
    do
        echo "calculating_omega_by_yn00.pl --omega_for_PSG 1.0 --target_branch_species $y ../orthologGroups_CDS/$i.fasta.align ../orthologGroups_CDS/$i.tree > ../dnds_yn00/$i.txt 2> dnds_yn00/$i.stats"
    done > command.calculating_omega_by_yn00.list
    ParaFly -c command.calculating_omega_by_yn00.list -CPU 20

    mkdir dnds_ML
    for i in `cat ../orthologGroups.txt`
    do
        echo "calculating_omega_by_codeml.pl --omega_for_PSG 1.0 --target_branch_species $y ../orthologGroups_CDS/$i.phy ../orthologGroups_CDS/$i.tree > ../dnds_ML/$i.txt 2> dnds_ML/$i.stats"
    done > command.calculating_omega_by_ML.list
    ParaFly -c command.calculating_omega_by_ML.list -CPU 20

    perl -e 'while (<dnds*/*.stats>) { my $cog = $1 if m/(OCG\d+)/; open IN, $_; <IN>; $_ = <IN>; @_ = split /\s+/; print "$cog\n" if $_[1] >= 1; close IN; }' | sort | uniq > candidate_PSG.list
    cd ..
done
# 统计各分枝的候选正选择基因数量
for x in `cat branches_outerNode.txt branches_innerNode.txt`
do
    x=${x/:*/}
    wc -l BS_$x/candidate_PSG.list | perl -e 'while (<>) { print "$2\t$1\n" if m/(\d+)\s+BS_(\w+)/; }'
done > results.candidate_PSG_number
# (3) 使用branch-site model对候选基因进行正选择分析
# 对内部节点分枝进行进行正选择基因分析
for x in `cat branches_innerNode.txt`
do
    y=$x
    x=${x/:*/}
    y=${y/*:/}

    mkdir BS_$x/branch-site_model
    for i in `cat BS_$x/candidate_PSG.list`
    do
        echo "paml_branch-site_model_analysis.pl --target_branch_species $y orthologGroups_CDS/$i.phy orthologGroups_CDS/$i.tree BS_$x/branch-site_model/$i > BS_$x/branch-site_model/$i.p"
    done >> command.paml_branch-site_model.list
done
ParaFly -c command.paml_branch-site_model.list -CPU 8
# 统计各分枝的正选择基因信息
for x in `cat branches_innerNode.txt`
do
    y=$x
    x=${x/:*/}
    y=${y/*:/}
    cd BS_$x/
    perl -e 'while (<branch-site_model/*.p>) { my $ocg = $_; $ocg =~ s/.*\///; $ocg =~ s/.p//; open IN, $_; my $p = "NULL\n"; while (<IN>) { chomp; @_ = split /\t/; $p = "$_[0]\n" if ($_[0] <= 0.05 && ($_[1] eq "BEB_Significance_YES")); } print "$ARGV[0]\t$ocg\t$p"; }' $x
    cd ..
done > results.PSG
# 对特定的外部节点分枝进行正选择基因分析
cd BS_laame; mkdir branch-site_model
for i in `cat candidate_PSG.list`
do
    echo "paml_branch-site_model_analysis.pl --target_branch_species sccit --adjacent_species_of_single_outNode parub ../orthologGroups_CDS/$i.phy ../orthologGroups_CDS/$i.tree branch-site_model/$i > branch-site_model/$i.p"
done > command.paml_branch-site_model.list
ParaFly -c command.paml_branch-site_model.list -CPU 5
perl -e 'while (<branch-site_model/*.p>) { my $ocg = $_; $ocg =~ s/.*\///; $ocg =~ s/.p//; open IN, $_; my $p = "NULL\n"; while (<IN>) { chomp; @_ = split /\t/; $p = "$_[0]\n" if ($_[0] <= 0.05 && ($_[1] eq "BEB_Significance_YES")); } print "$ARGV[0]\t$ocg\t$p"; }' sccit > results.PSG
cd ..

cd ..


## g. 使用MCScanX进行共线性区块分析
mkdir -p /home/train/14.genome_comparison/g.MCScanX
cd /home/train/14.genome_comparison/g.MCScanX

# 准备2个物种基因组的蛋白质序列文件和GFF文件
ln -s ../a.preparing_data/laame.geneModels.gff3 ./
ln -s ../a.preparing_data/laame.pep.fasta ./
ln -s ../a.preparing_data/laame.genome.fasta ./
ln -s ../a.preparing_data/plost.geneModels.gff3 .
ln -s ../a.preparing_data/plost.pep.fasta .
ln -s ../a.preparing_data/plost.genome.fasta ./

cat laame.pep.fasta plost.pep.fasta > all.fasta
perl -p -i -e 's/\*$//; s/\*/X/g;' all.fasta 
diamond makedb --in all.fasta --db all
diamond blastp --db all --query all.fasta --out diamond.out --outfmt 5 --sensitive --max-target-seqs 100 --evalue 1e-5 --id 10 --tmpdir /dev/shm --threads 8
mkdir data
parsing_blast_result.pl --no-header --max-hit-num 100 --evalue 1e-6 --identity 0.5 --subject-coverage 0.5 --query-coverage 0.5 diamond.out > data/input.blast
perl -e 'while (<>) { if (m/^(\S+)\t.*\tgene\t(\d+)\t(\d+).*ID=([^\s;]+)/) { print "$1\t$4\t$2\t$3\n" } }' plost.geneModels.gff3 laame.geneModels.gff3 > data/input.gff
MCScanX data/input
grep -v -P "plost.*plost" data/input.collinearity | grep -P "plost" > data/input.collinearity_interspecific
mcscanx_stats_blocks.pl data/input.collinearity_interspecific > data/input.collinearity_interspecific.stats
grep -P "plost.*plost" data/input.collinearity > data/input.collinearity_intraspecific
mcscanx_stats_blocks.pl data/input.collinearity_intraspecific > data/input.collinearity_intraspecific.stats

circos_from_MCScanX_out.pl --out-ref-WGD out_WGD --ref-gff3 plost.geneModels.gff3 --ref-fasta plost.genome.fasta --query-gff3 laame.geneModels.gff3 --query-fasta laame.genome.fasta --ref-label PO --query-label LA --min-block-size 5 data/input.collinearity
cd out
circos -conf circos.conf
cd ../out_WGD
ircos -conf circos.conf
cd ..


## h. 使用Mummer对两个基因组进行比较
mkdir -p /home/train/14.genome_comparison/Mummer
cd /home/train/14.genome_comparison/Mummer

cp ~/00.incipient_data/data_for_genome_assembling/assemblies_of_Malassezia_sympodialis/Malassezia_sympodialis.genome_V01.fasta ./
cp ~/00.incipient_data/data_for_genome_assembling/assemblies_of_Malassezia_sympodialis/IDBA.fasta ./
cp ~/00.incipient_data/data_for_genome_assembling/assemblies_of_Malassezia_sympodialis/MaSuRCA.fasta .

para_nucmer --CPU 8 Malassezia_sympodialis.genome_V01.fasta IDBA.fasta > out.delta
#nucmer -c 200 -g 200 -p out Malassezia_sympodialis.genome_V01.fasta IDBA.fasta
delta-filter -i 95 -r -q out.delta > out.rq.delta
show-coords -c -d -l -I 95 -L 10000 -r out.rq.delta > out.show
mummerplot -f -l -p out -s large -t png out.delta
gnuplot out.gp
