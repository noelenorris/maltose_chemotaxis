clear all;
clc;

all_results = [];

file_id = fopen('results.txt')

result = fgetl(file_id)

while(result ~= -1)
    load(result)
    all_results = [all_results; list_results];
    result = fgetl(file_id)
end

fclose(file_id);
