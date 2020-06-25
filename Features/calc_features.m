function output = calc_features(RNA,ORF_start,UTR3_start,region,seed_length,RNA_start,miRNA,seed_type,...
    phastcons100,phastcons20,phylops100,phylops20)
%Calculate a 1xn table of features

addpath('../Utils')

UTR5 = RNA(1:ORF_start-1);
ORF = RNA(ORF_start:UTR3_start-1);
UTR3 = RNA(UTR3_start:end);

miRNA = upper(miRNA);
miRNA(miRNA == 'T') = 'U';
ORF_start = length(UTR5) + 1;

if any(strcmp(seed_type,{'7mer-A1','8mer'}))
    RNA_pos1 = RNA_start + seed_length-1;
else
    RNA_pos1 = RNA_start + seed_length;
end

%region_start/end - in local (RNA strand) coordinates
switch region
    case 'UTR5'
        region_start = 1;
        region_end = length(UTR5);
    case 'ORF'
        region_start = length(UTR5)+1;
        region_end = length(UTR5)+length(ORF);
    case 'UTR3'
        region_start = length(UTR5)+length(ORF)+1;
        region_end = length(RNA);
end

output = table(0);

%Conservation
cd 'Conservation';

output.prob_binom = calc_prob_binom(RNA,RNA_start,seed_length);
output.prob_exact = calc_prob_exact(RNA,RNA_start,seed_length);

[output.phastcons_seed100,output.phastcons_site100,output.phastcons_flank100,phylops_seed100,output.phylops_site100...
    ,output.phylops_flank100,output.phastcons_seed20,output.phastcons_site20,output.phastcons_flank20,phylops_seed20,...
    output.phylops_site20,output.phylops_flank20] = calc_cons(phastcons100,phastcons20,phylops100,phylops20,RNA_start,seed_length);

phylops_seed100 = flipud(phylops_seed100);
for g = 1:length(phylops_seed100)
    output{1,sprintf('phylops100_pos%d',g)} = phylops_seed100(g);
end

phylops_seed20 = flipud(phylops_seed20);
for g = 1:length(phylops_seed20)
    output{1,sprintf('phylops20_pos%d',g)} = phylops_seed20(g);
end

%Sequence
cd '../Sequence';

[output.distance_score, output.dist_start, output.dist_end] = calc_utr_pos(RNA_start,region_start,region_end,miRNA,seed_length);
output.au_content = calc_au_content(RNA,RNA_start,seed_length,seed_type);
output.miR_pos1_A = double(miRNA(1) == 'A');
output.miR_pos1_C = double(miRNA(1) == 'C');
output.miR_pos1_G = double(miRNA(1) == 'G');
output.miR_pos8_A = double(miRNA(8) == 'A');
output.miR_pos8_C = double(miRNA(8) == 'C');
output.miR_pos8_G = double(miRNA(8) == 'G');

if RNA_pos1 <= length(RNA)
    output.RNA_pos1_A = double(RNA(RNA_pos1) == 'A');
    output.RNA_pos1_C = double(RNA(RNA_pos1) == 'C');
    output.RNA_pos1_G = double(RNA(RNA_pos1) == 'G');
end

RNA_pos8 = RNA_pos1 - 7;
if RNA_pos8 >= 1
    output.RNA_pos8_A = double(RNA(RNA_pos8) == 'A');
    output.RNA_pos8_C = double(RNA(RNA_pos8) == 'C');
    output.RNA_pos8_G = double(RNA(RNA_pos8) == 'G');
    if RNA_pos8 >= 2
        output.RNA_pos9_A = double(RNA(RNA_pos8-1) == 'A');
        output.RNA_pos9_C = double(RNA(RNA_pos8-1) == 'C');
        output.RNA_pos9_G = double(RNA(RNA_pos8-1) == 'G');
    else
        output.RNA_pos9_A = 0;
        output.RNA_pos9_C = 0;
        output.RNA_pos9_G = 0;
    end
else
    output.RNA_pos8_A = 0;
    output.RNA_pos8_C = 0;
    output.RNA_pos8_G = 0;
    output.RNA_pos9_A = 0;
    output.RNA_pos9_C = 0;
    output.RNA_pos9_G = 0;
end

switch region
    case 'UTR5'
        curr_region = UTR5;
    case 'ORF'
        curr_region = ORF;
    case 'UTR3'
        curr_region = UTR3;
end

output.dist_start_rel = round(10 ^ output.dist_start) / length(curr_region);
output.region_A = count(curr_region,'A') / length(curr_region);
output.region_C = count(curr_region,'C') / length(curr_region);
output.region_G = count(curr_region,'G') / length(curr_region);

seed_seq = RNA(RNA_start:RNA_start + seed_length-1);
output.site_A = count(seed_seq,'A') / seed_length;
output.site_C = count(seed_seq,'C') / seed_length;
output.site_G = count(seed_seq,'G') / seed_length;

if ~isempty(UTR5)
    output.UTR5_len = log10(length(UTR5));
    output.UTR5_AU = 1 - count(UTR5,{'G','C'}) / length(UTR5);
else
    output.UTR5_len = 0;
    output.UTR5_AU = 0;
end

output.ORF_len = log10(length(ORF));
output.ORF_AU = 1 - count(ORF,{'G','C'}) / length(ORF);

if ~isempty(UTR3)
    output.UTR3_len = log10(length(UTR3));
    output.UTR3_AU = 1 - count(UTR3,{'G','C'}) / length(UTR3);
else
    output.UTR3_len = 0;
    output.UTR3_AU = 0;
end

load('../../Data/TS7.mat','TS7');
ind = find(strcmp(TS7.seed,miRNA(2:8)));
output.TA = TS7.TA(ind);
output.TA_hela = TS7.TA_hela(ind);
output.SPS_6mer = TS7.SPS_6mer(ind);
output.SPS_7mer_m8 = TS7.SPS_7mer_m8(ind);
output.mean_SPS = TS7.mean_SPS(ind);

miR_end = max(1,RNA_pos1 - (length(miRNA) - 1));
pairing = erase(num2str(is_watson_crick(RNA(miR_end:miR_end+length(miRNA)-1),miRNA)),' ');
tokens = regexp(pairing,'1+','match');
if ~isempty(tokens)
    output.longest_pairing = max(cellfun(@length,tokens));
else
    output.longest_pairing = 0;
end

output.paired_miR_end = sum(is_watson_crick(RNA(miR_end:miR_end+6),miRNA(end-6:end)));
output.pairing3p_sw = calc_3p_pairing_sw(RNA,RNA_start,miRNA,seed_type);
output.pairing3p = calc_3p_pairing(RNA,RNA_start,miRNA,seed_length,2);

miR_cog = RNA(miR_end:RNA_pos1);
output.miR_cog_A = count(miR_cog,'A') / length(miR_cog);
output.miR_cog_C = count(miR_cog,'C') / length(miR_cog);
output.miR_cog_G = count(miR_cog,'G') / length(miR_cog);

miR_cog_up = RNA(max(1,miR_end-50):max(1,miR_end-1));
output.miR_cog_up_A = count(miR_cog_up,'A') / length(miR_cog_up);
output.miR_cog_up_C = count(miR_cog_up,'C') / length(miR_cog_up);
output.miR_cog_up_G = count(miR_cog_up,'G') / length(miR_cog_up);
miR_cog_down = RNA(min(length(RNA),RNA_pos1+1):min(length(RNA),RNA_pos1+50));
output.miR_cog_down_A = count(miR_cog_down,'A') / length(miR_cog_down);
output.miR_cog_down_C = count(miR_cog_down,'C') / length(miR_cog_down);
output.miR_cog_down_G = count(miR_cog_down,'G') / length(miR_cog_down);

[ORF_pum,UTR3_pum,nei_pum] = calc_pum(RNA,RNA_start,seed_length,ORF_start,UTR3_start-1);
output.ORF_pum1 = ORF_pum(1);
output.ORF_pum2 = ORF_pum(2);
output.UTR3_pum1 = UTR3_pum(1);
output.UTR3_pum2 = UTR3_pum(2);
output.nei_pum1 = nei_pum(1);
output.nei_pum2 = nei_pum(2);

RBP_sites = {'UAUUU','UAUUUAU','AUUUAU','GUGAAG','AUUUAUU','UAUUUA','UAUUUU','AGCCA','AGAGAA','UUUAU','UUUUU','UUUUAAA','UUUUUU',...
'UUUUUUU','UUUAAA','UUUAAAA','UUUUG','UGUUU','GAAAA','CAAAA'};

for j = 1:length(RBP_sites)
    if ~isempty(UTR3)
        output{1,sprintf('RBP%d',j)} = count(UTR3,RBP_sites{j});
    else
        output{1,sprintf('RBP%d',j)} = 0;
    end
end

seed_table = find_potential_targets(RNA,miRNA,ORF_start,UTR3_start,0);

for seed = ["6mer-A1","offset-7mer","offset-6mer","centered","CDNST1","CDNST2","CDNST3","CDNST4"]
    output{1,sprintf('utr3_%s',strrep(seed,'-','_'))} = sum(strcmp(seed_table.region,'UTR3') & strcmp(seed_table.seed_type,seed));
end

miR_T = seqrcomplement(miRNA(2:4));
miR_T(miR_T == 'U') = 'T';

if ~isempty(UTR5)
    UTR5_T = UTR5;
    UTR5_T(UTR5_T == 'U') = 'T';
    output.utr5_3mer = length(regexp(UTR5_T,miR_T));
else
    output.utr5_3mer = 0;
end

if ~isempty(ORF)
    ORF_T = ORF;
    ORF_T(ORF_T == 'U') = 'T';
    output.orf_3mer = length(regexp(ORF_T,miR_T));
else
    output.orf_3mer = 0;
end

if ~isempty(UTR3)
    UTR3_T = UTR3;
    UTR3_T(UTR3_T == 'U') = 'T';
    output.utr3_3mer = length(regexp(UTR3_T,miR_T));
else
    output.utr3_3mer = 0;
end

%Thermodynamic
cd '../Thermo';

[output.dg_duplex,output.dg_duplex_seed,output.dg_binding,output.dg_binding_seed] = ...
    calc_dg_duplex_binding(RNA,miRNA,2,RNA_pos1,RNA_start,seed_type,seed_length); %#ok<*AGROW>

output.dg_open = calc_dg_open(RNA,RNA_pos1,length(miRNA));

output.SPS = calc_sps(RNA,miRNA,RNA_start,seed_type,seed_length);

sa = calc_sa(RNA,RNA_start,seed_type,region_end,region);
if ~isnan(sa)
    output.SA = sa;
else
    output.SA = 999;
end

max_energy = 0:-5:-30;
mismatch = 2;

MREs = calc_MRE(ORF,miRNA,max_energy,mismatch);
for p = 1:length(MREs)
    output{1,sprintf('MRE_therm_ORF_m2_%d',abs(max_energy(p)))} = MREs(p);
end

MREs = calc_MRE(UTR3,miRNA,max_energy,mismatch);
for p = 1:length(MREs)
    output{1,sprintf('MRE_therm_UTR_m2_%d',abs(max_energy(p)))} = MREs(p);
end

% Translation
cd '../Translation';

[output.cub_global_win,output.cub_global_orf,output.cub_local_win,output.cub_local_orf,output.cub_local_CAI_win...
,output.cub_local_CAI_orf,output.tAI_score_win,output.tAI_score_orf,output.charge_score_win...
,output.charge_score_orf,output.slow_score_win,output.slow_score_orf,output.TDR_score_win,output.TDR_score_orf] = ...
calc_translation(RNA,RNA_start,seed_length,ORF_start,UTR3_start-1,region);

cd '..'

output.Var1 = [];

end