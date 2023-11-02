function [new_datax, new_datay, new_s_datay] = avoidOversampling(datax, datay, sigmay)

    
    % safety check
    assert(length(datax) == length(datay), "datax and datay should have the same length but length(datax) = '" + length(datax) + "' and length(datay) = '" + length(datay) + "'");

    % vogliamo che i dati in X siano ordinati in modo crescente
    [datax, order] = sort(datax);
    datay = datay(order);

    new_datax = [];
    new_datay = [];
    new_s_datay = [];
    last_x = datax(1);
    start_cur_list = 1;
    end_cur_list = 1;
    for ii = 1:length(datax)
        assert(datax(ii) >= last_x, "sorted_datax should be a non-decreasing sequence.");
        if datax(ii) == last_x
            end_cur_list = ii;
        else
            if(start_cur_list == end_cur_list)
                new_datax(length(new_datax) + 1, 1) = last_x;
                new_datay(length(new_datay) + 1, 1) = datay(start_cur_list);
                new_s_datay(length(new_s_datay) + 1, 1) = sigmay(start_cur_list);  
            else    
                tmp_sigma = std(datay(start_cur_list:end_cur_list));
                if(tmp_sigma == 0)
                    new_datax(length(new_datax) + 1:length(new_datax) + 1 + end_cur_list - start_cur_list, 1) = datax(start_cur_list:end_cur_list);
                    new_datay(length(new_datay) + 1:length(new_datay) + 1 + end_cur_list - start_cur_list, 1) = datay(start_cur_list:end_cur_list);
                    new_s_datay(length(new_s_datay) + 1:length(new_s_datay) + 1 + end_cur_list - start_cur_list, 1) = sigmay(start_cur_list:end_cur_list);
                else
                    new_datax(length(new_datax) + 1, 1) = last_x;
                    new_datay(length(new_datay) + 1, 1) = mean(datay(start_cur_list:end_cur_list));
                    new_s_datay(length(new_s_datay) + 1, 1) = tmp_sigma;    
                end
            end
            start_cur_list = ii;
            end_cur_list = ii;
            last_x = datax(ii);
        end
    end
    if start_cur_list ~= end_cur_list
        tmp_sigma = std(datay(start_cur_list:end_cur_list));
        if(tmp_sigma == 0)
            new_datax(length(new_datax) + 1:length(new_datax) + 1 + end_cur_list - start_cur_list, 1) = datax(start_cur_list:end_cur_list);
            new_datay(length(new_datay) + 1:length(new_datay) + 1 + end_cur_list - start_cur_list, 1) = datay(start_cur_list:end_cur_list);
            new_s_datay(length(new_s_datay) + 1:length(new_s_datay) + 1 + end_cur_list - start_cur_list, 1) = sigmay(start_cur_list:end_cur_list);
        else
            new_datax(length(new_datax) + 1, 1) = last_x;
            new_datay(length(new_datay) + 1, 1) = mean(datay(start_cur_list:end_cur_list));
            new_s_datay(length(new_s_datay) + 1, 1) = tmp_sigma;    
        end
    end
end