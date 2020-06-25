function seed_table = find_potential_targets(RNA,miRNA,ORF_start,UTR3_start,only_canonical)

%Find potential miRNA targets in RNA

%RNA = The sequence of the mRNA, including UTRs [char]

%miRNA = The sequence of the miRNA [char]

%ORF/UTR3_start = The first coordinate of the ORF/UTR3 region in "mRNA"
%[double]

%only_canonical = 0/1 index to find only canonical target sites, or both
%canonical and non-canonical sites [double]

RNA = upper(RNA);
RNA(RNA == 'U') = 'T';
miRNA = upper(miRNA);
miRNA(miRNA == 'U') = 'T';

seed_table = table(0);

sites_8mer = regexp(RNA,[seqrcomplement(miRNA(8)),'(?=',seqrcomplement(miRNA(2:7)),'A)'])';
if ~isempty(sites_8mer)
    temp_table = table(sites_8mer,cellstr(repmat('8mer',length(sites_8mer),1)),repmat(8,length(sites_8mer),1),...
        'VariableNames',{'RNA_start','seed_type','seed_length'});
    seed_table = temp_table;
end

sites_7mer_m8 = regexp(RNA,[seqrcomplement(miRNA(8)),'(?=',seqrcomplement(miRNA(2:7)),'[^A])'])';
if ~isempty(sites_7mer_m8)
    temp_table = table(sites_7mer_m8,cellstr(repmat('7mer-m8',length(sites_7mer_m8),1)),repmat(7,length(sites_7mer_m8),1),...
        'VariableNames',{'RNA_start','seed_type','seed_length'});
    if width(seed_table) == 1
        seed_table = temp_table;
    else
        seed_table = [seed_table; temp_table];
    end
end

sites_7mer_A1 = regexp(RNA,[sprintf('[^%s](?=',seqrcomplement(miRNA(8))),seqrcomplement(miRNA(2:7)),'A)'])';
if ~isempty(sites_7mer_A1)
    temp_table = table(sites_7mer_A1+1,cellstr(repmat('7mer-A1',length(sites_7mer_A1),1)),repmat(7,length(sites_7mer_A1),1),...
        'VariableNames',{'RNA_start','seed_type','seed_length'});
    if width(seed_table) == 1
        seed_table = temp_table;
    else
        seed_table = [seed_table; temp_table];
    end
end

sites_6mer = regexp(RNA,[sprintf('[^%s](?=',seqrcomplement(miRNA(8))),seqrcomplement(miRNA(2:7)),'[^A])'])';
if ~isempty(sites_6mer)
    temp_table = table(sites_6mer+1,cellstr(repmat('6mer',length(sites_6mer),1)),repmat(6,length(sites_6mer),1),...
        'VariableNames',{'RNA_start','seed_type','seed_length'});
    if width(seed_table) == 1
        seed_table = temp_table;
    else
        seed_table = [seed_table; temp_table];
    end
end

if ~only_canonical
    sites_6mer_A1 = regexp(RNA,[sprintf('[^%s](?=',seqrcomplement(miRNA(7))),seqrcomplement(miRNA(2:6)),'A)'])';
    if ~isempty(sites_6mer_A1)
        temp_table = table(sites_6mer_A1+1,cellstr(repmat('6mer-A1',length(sites_6mer_A1),1)),repmat(5,length(sites_6mer_A1),1),...
            'VariableNames',{'RNA_start','seed_type','seed_length'});
        if width(seed_table) == 1
            seed_table = temp_table;
        else
            seed_table = [seed_table; temp_table];
        end
    end

    sites_offset_7mer = regexp(RNA,[seqrcomplement(miRNA(9)),'(?=',...
        seqrcomplement(miRNA(3:8)),sprintf('[^%s]',seqrcomplement(miRNA(2))),')'])';
    if ~isempty(sites_offset_7mer)
        temp_table = table(sites_offset_7mer,cellstr(repmat('offset-7mer',length(sites_offset_7mer),1)),...
            repmat(7,length(sites_offset_7mer),1),'VariableNames',{'RNA_start','seed_type','seed_length'});
        if width(seed_table) == 1
            seed_table = temp_table; %#ok<*NASGU>
        else
            seed_table = [seed_table; temp_table];
        end
    end


    sites_offset_6mer = regexp(RNA,[sprintf('[^%s](?=',...
        seqrcomplement(miRNA(9))),seqrcomplement(miRNA(3:8)),sprintf('[^%s]',seqrcomplement(miRNA(2))),')'])';
    if ~isempty(sites_offset_6mer)
        temp_table = table(sites_offset_6mer+1,cellstr(repmat('offset-6mer',length(sites_offset_6mer),1)),...
            repmat(6,length(sites_offset_6mer),1),'VariableNames',{'RNA_start','seed_type','seed_length'});
        if width(seed_table) == 1
            seed_table = temp_table;
        else
            seed_table = [seed_table; temp_table];
        end
    end

    sites_centered = regexp(RNA,[seqrcomplement(miRNA(13)),'(?=',seqrcomplement(miRNA(3:12)),')|',seqrcomplement(miRNA(14)),'(?=',seqrcomplement(miRNA(4:13)),')|',seqrcomplement(miRNA(15)),'(?=',seqrcomplement(miRNA(5:14)),')'])';
    if ~isempty(sites_centered)
        temp_table = table(sites_centered+1,cellstr(repmat('centered',length(sites_centered),1)),...
            repmat(5,length(sites_centered),1),'VariableNames',{'RNA_start','seed_type','seed_length'});
        if width(seed_table) == 1
            seed_table = temp_table;
        else
            seed_table = [seed_table; temp_table];
        end
    end
    
    sites_CDNST1 = regexp(RNA,[sprintf('[^%s](?=',seqrcomplement(miRNA(7))),seqrcomplement(miRNA(2:6)),'[^A])'])';
    if ~isempty(sites_CDNST1)
        temp_table = table(sites_CDNST1+1,cellstr(repmat('CDNST1',length(sites_CDNST1),1)),...
            repmat(5,length(sites_CDNST1),1),'VariableNames',{'RNA_start','seed_type','seed_length'});
        if width(seed_table) == 1
            seed_table = temp_table;
        else
            seed_table = [seed_table; temp_table];
        end
    end

    sites_CDNST2 = regexp(RNA,[seqrcomplement(miRNA(7)),'(?=',...
        seqrcomplement(miRNA(6)),sprintf('[^%s]',seqrcomplement(miRNA(5))),seqrcomplement(miRNA(2:4)),'A)'])';
    if ~isempty(sites_CDNST2)
        temp_table = table(sites_CDNST2,cellstr(repmat('CDNST2',length(sites_CDNST2),1)),...
            repmat(6,length(sites_CDNST2),1),'VariableNames',{'RNA_start','seed_type','seed_length'});
        if width(seed_table) == 1
            seed_table = temp_table;
        else
            seed_table = [seed_table; temp_table];
        end
    end

    sites_CDNST3 = regexp(RNA,[seqrcomplement(miRNA(9)),'(?=',seqrcomplement(miRNA(7:8)),sprintf('[^%s]',seqrcomplement(miRNA(6))),seqrcomplement(miRNA(5)),...
        sprintf('[^%s]',seqrcomplement(miRNA(4))),seqrcomplement(miRNA(3)),sprintf('[^%s]',seqrcomplement(miRNA(2))),')'])';
    if ~isempty(sites_CDNST3)
        temp_table = table(sites_CDNST3,cellstr(repmat('CDNST3',length(sites_CDNST3),1)),...
            repmat(7,length(sites_CDNST3),1),'VariableNames',{'RNA_start','seed_type','seed_length'});
        if width(seed_table) == 1
            seed_table = temp_table;
        else
            seed_table = [seed_table; temp_table];
        end
    end

    sites_CDNST4 = regexp(RNA,[seqrcomplement(miRNA(8)),'(?=',sprintf('[^%s]',seqrcomplement(miRNA(7))),...
        sprintf('[^%s]',seqrcomplement(miRNA(6))),sprintf('[^%s]',seqrcomplement(miRNA(5))),seqrcomplement(miRNA(2:4)),')'])';
    if ~isempty(sites_CDNST4)
        temp_table = table(sites_CDNST4,cellstr(repmat('CDNST4',length(sites_CDNST4),1)),...
            repmat(3,length(sites_CDNST4),1),'VariableNames',{'RNA_start','seed_type','seed_length'});
        if width(seed_table) == 1
            seed_table = temp_table;
        else
            seed_table = [seed_table; temp_table];
        end
    end
end

if width(seed_table) == 1
    seed_table(1,:) = [];
else
    seed_table{seed_table.RNA_start < ORF_start,'region'} = {'UTR5'};
    seed_table{seed_table.RNA_start >= ORF_start & seed_table.RNA_start < UTR3_start ,'region'} = {'ORF'};
    seed_table{seed_table.RNA_start >= UTR3_start,'region'} = {'UTR3'};
end

end