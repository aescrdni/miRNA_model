function output = predict(RNA,ORF_start,UTR3_start,miRNA,RNA_start,...
    seed_type,phastcons20,phastcons100,phylops20,phylops100,riboseq,with_features)

%Calculating features for a binding site table
%RNA = The sequence of the mRNA, including UTRs [char]

%ORF/UTR3_start = The first coordinate of the ORF/UTR3 region in "mRNA"
%[double]

%miRNA = The sequence of the miRNA [char]

%RNA_start = The first coordinate of the target site in the mRNA (i.e. the
%one paired with position 7 in the miRNA for 6mer/7mer-A1 sites, or the one
%paired with position 8 in the miRNA for 7mer-m8/8mer sites) [double]

%seed_type = 6mer/7mer-A1/7mer-m8/8mer [char]

%phastcons20/100 - A vector of the PhastCons scores for each coordinate of
%the mRNA, based on multiple alignment of hg38 with 99 vertebrates/19
%mammals, respectively (can be downloaded from the UCSC Genome Browser's database)
%[double vector of size 1x(mRNA length)]

%phylops20/100 - Similar to phastcons20/100, based on PhyloP scores [double vector of size 1x(mRNA length)]

%riboseq - A 1x3 vector, containing average RiboSeq score for the first,
%second and last third of the ORF [1x3 double vector]

%with_features - 0/1 index for including calculated features in the output or only the
%predicted repression [double]

seeds = {'6mer','7mer-A1','7mer-m8','8mer'};
if isempty(find(strcmp(seeds,seed_type),1))
    error('seed_type should be 6mer/7mer-A1/7mer-m8/8mer');
end

seed_length = str2double(seed_type(1));

cd 'Features'

RNA = upper(RNA);
RNA(RNA == 'T') = 'U';

miRNA = upper(miRNA);
miRNA(miRNA == 'T') = 'U';

if RNA_start < ORF_start
    region = 'UTR5';
elseif RNA_start > UTR3_start
    region = 'UTR3';
else
    region = 'ORF';
end

%Calculating features
output = calc_features(RNA,ORF_start,UTR3_start,region,seed_length,RNA_start,miRNA,seed_type,...
    phastcons100,phastcons20,phylops100,phylops20);

for t = 1:3
    output{1,sprintf('RiboSeq_third%d',t)} = riboseq(t);
end

cd '..'

%Loading model and predicting repression

load('Data/models.mat','models');

if strcmp(region,'UTR3')
    r = 2;
else
    r = 1;
end

s = find(strcmp(seeds,seed_type));

model = models{r,s};

features_mat = table2array(output(1,model.features));
repress = model.intercept + features_mat * model.coeffs;

if r == 1 || s == 1
    repress = min(repress,0);
elseif s == 2
    repress = min(repress,-0.01);
elseif s == 3
    repress = min(repress,-0.02);
else
    repress = min(repress,-0.03);
end

if with_features
    output.repression(1) = repress;
else
    output = repress;
end

end