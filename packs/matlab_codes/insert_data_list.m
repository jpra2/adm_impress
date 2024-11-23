function [data_list, count_data] = insert_data_list(old_data_list, old_count, data_to_insert)

data_list = old_data_list;
data_list(old_count, :) = data_to_insert;
count_data = old_count + 1;

end