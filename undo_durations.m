reprocess_data

for Mind = 1:1
    for Pind = 1:5
        for para = 1:1
            for i = 1:9
                t= test_data_tau5_noeps39_lambda0001{Pind,i};
                test_u_end = [u_total(i,:),0];
                u_star_orig = t.full_deconvolved_BrAC;
                u_star_end = u_star_orig;
                [peak_est, peaktime_est] = max(u_star_end);
                [peak_act, peaktime_act] = max(test_u_end);


                t.actual_error=u_star_end-test_u_end;
                t.L2_error=sum((u_star_end-test_u_end).^2);
                t.Linf_error=max(abs(u_star_end-test_u_end));
                t.AUC_error=sum(u_star_end)-sum(test_u_end);
                t.peak_time_error=tau*(peaktime_est-peaktime_act);
                t.peak_height_error=peak_est-peak_act;

                test_data_tau5_noeps39_lambda0001{Mind,Pind,para,i}=t;
            end
        end
    end
end
