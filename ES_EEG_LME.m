%% LME analysis for first two bins 

cd D:\ES_data_Hannah\for_stats

temp = dir('*'); % assuming 3 directories

for i=1:length(temp)

    if temp(i).isdir == 1 & length(temp(i).name) > 3
        cd(temp(i).name)

        load for_stats.mat
        load for_stats_expquarter.mat
        load for_stats_subject.mat
        load for_stats_validinvalid.mat

       
        % do ANOVA as usual

        [p,t,stats,terms] = anovan(for_stats(:),{for_stats_expquarter(:) for_stats_validinvalid(:) for_stats_subject(:)}, ...
            'varnames',{'qt' 'vi' 'sb'}, 'random',3,'model','interaction','display','off');

        display(' ')
        display('*****')
        display(temp(i).name)
        display('*****')
        display(' ')
        display('RM ANOVA, all data points')
        display(strcat(['Main effect of quarter, F(1,29) = ',num2str(cell2mat(t(2,6))), ', p = ',num2str(cell2mat(t(2,7)))]))
        display(strcat(['Main effect of validity, F(1,29) = ',num2str(cell2mat(t(3,6))), ', p = ',num2str(cell2mat(t(3,7)))]))
        display(strcat(['Interaction quarter x validity, F(1,29) = ',num2str(cell2mat(t(5,6))), ', p = ',num2str(cell2mat(t(5,7)))]))
        display(' ')

        % fit linear mixed model, remove outliers, and re-fit

        tbl = table(for_stats(:),for_stats_expquarter(:),for_stats_validinvalid(:),categorical(for_stats_subject(:)), ...
            'VariableNames',{'acc' 'qtr' 'val' 'sub'});

        mdl = fitlm(tbl,'acc ~ 1 + qtr + val + qtr*val + sub');

        display('LME, all data points')
        display(strcat(['Main effect of quarter, T(',num2str(mdl.DFE),') = ',num2str(table2array(mdl.Coefficients(2,3))), ', p = ',num2str(table2array(mdl.Coefficients(2,end)))]))
        display(strcat(['Main effect of validity, T(',num2str(mdl.DFE),') = ',num2str(table2array(mdl.Coefficients(3,3))), ', p = ',num2str(table2array(mdl.Coefficients(3,end)))]))
        display(strcat(['Interaction quarter x validity, T(',num2str(mdl.DFE),') = ',num2str(table2array(mdl.Coefficients(end,3))), ', p = ',num2str(table2array(mdl.Coefficients(end,end)))]))
        display(' ')


        outliers = find((mdl.Diagnostics.CooksDistance)>3*mean(mdl.Diagnostics.CooksDistance));
        tbl(outliers,:)=[];
        mdl2 = fitlm(tbl,'acc ~ 1 + qtr + val + qtr*val + sub');

        R2_model = mdl2.Rsquared.Ordinary;
        cohensf = sqrt(R2_model / (1 - R2_model));


        display('LME, no outliers')
        display(strcat(['Main effect of quarter, T(',num2str(mdl2.DFE),') = ',num2str(table2array(mdl2.Coefficients(2,3))), ', p = ',num2str(table2array(mdl2.Coefficients(2,end)))]))
        display(strcat(['Main effect of validity, T(',num2str(mdl2.DFE),') = ',num2str(table2array(mdl2.Coefficients(3,3))), ', p = ',num2str(table2array(mdl2.Coefficients(3,end)))]))
        display(strcat(['Interaction quarter x validity, T(',num2str(mdl2.DFE),') = ',num2str(table2array(mdl2.Coefficients(end,3))), ', p = ',num2str(table2array(mdl2.Coefficients(end,end)))]))
        display(strcat(['Cohen''s f for the LME model (no outliers) = ', num2str(cohensf)]));
        display(' ')
        cd ..

        if table2array(mdl2.Coefficients(end,end))<.05
            tbl1 = tbl(find(tbl.qtr==1),[1 3 4]);
            mdl_qt1 = fitlm(tbl1,'acc ~ 1 + val + sub');

            mdl_qt1_no_val = fitlm(tbl1, 'acc ~ 1 + sub'); 
            R2_qt1_full = mdl_qt1.Rsquared.Ordinary;
            R2_qt1_no_val = mdl_qt1_no_val.Rsquared.Ordinary;
            cohens_f_qt1_val = (R2_qt1_full - R2_qt1_no_val) / (1 - R2_qt1_full);

            display('LME, no outliers, post-hoc, QT1')
            display(strcat(['Main effect of validity, T(',num2str(mdl_qt1.DFE),') = ',num2str(table2array(mdl_qt1.Coefficients(2,3))), ', p = ',num2str(table2array(mdl_qt1.Coefficients(2,end)))]))
            display(strcat(['Partial Cohen''s f for validity, QT1 = ', num2str(cohens_f_qt1_val)]));
            display(' ')


            tbl2 = tbl(find(tbl.qtr==2),[1 3 4]);
            mdl_qt2 = fitlm(tbl2,'acc ~ 1 + val + sub');

            mdl_qt2_no_val = fitlm(tbl2, 'acc ~ 1 + sub'); 
            R2_qt2_full = mdl_qt2.Rsquared.Ordinary;
            R2_qt2_no_val = mdl_qt2_no_val.Rsquared.Ordinary;
            cohens_f_qt2_val = (R2_qt2_full - R2_qt2_no_val) / (1 - R2_qt2_full);

            display('LME, no outliers, post-hoc, QT2')
            display(strcat(['Main effect of validity, T(',num2str(mdl_qt2.DFE),') = ',num2str(table2array(mdl_qt2.Coefficients(2,3))), ', p = ',num2str(table2array(mdl_qt2.Coefficients(2,end)))]))
            display(strcat(['Partial Cohen''s f for validity, QT2 = ', num2str(cohens_f_qt2_val)]));
            display(' ')

        end


    end
end
