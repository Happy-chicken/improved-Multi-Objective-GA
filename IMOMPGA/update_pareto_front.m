function [pareto_front, pareto_set] = update_pareto_front(population, fitnesses, pareto_front, pareto_set)
% 更新帕累托前沿和帕累托集合

num_individuals = length(population);
num_pareto = size(pareto_front, 1);
new_individuals = [population, fitnesses];

% 将新个体加入帕累托集合
pareto_set = [pareto_set; population];
pareto_set_fitnesses = calc_fitnesses(pareto_set);

% 计算所有个体的适应度排名
[~, ranks] = sortrows(new_individuals(:,num_pareto+1:end));
fronts = cell(1,1);
fronts{1} = ranks(1:num_pareto);

% 生成帕累托前沿
for i = 1:num_pareto
    dominated_indices = [];
    for j = 1:num_individuals
        if j ~= i && dominates(new_individuals(j,num_pareto+1:end), new_individuals(i,num_pareto+1:end))
            dominated_indices = [dominated_indices, j];
        end
    end
    new_individuals(i,num_pareto+1) = length(dominated_indices);
    if new_individuals(i,num_pareto+1) == 0
        pareto_front = [pareto_front; new_individuals(i,num_pareto+2:end)];
    end
end

% 生成其他帕累托前沿
front_index = 1;
while ~isempty(fronts{front_index})
    next_front = [];
    for i = 1:length(fronts{front_index})
        dominated_indices = [];
        for j = 1:num_individuals
            if j ~= fronts{front_index}(i) && dominates(new_individuals(j,num_pareto+1:end), new_individuals(fronts{front_index}(i),num_pareto+1:end))
                dominated_indices = [dominated_indices, j];
            end
        end
        new_individuals(fronts{front_index}(i),num_pareto+1) = length(dominated_indices);
        if new_individuals(fronts{front_index}(i),num_pareto+1) == 0
            next_front = [next_front, fronts{front_index}(i)];
        end
    end
    front_index = front_index + 1;
    fronts{front_index} = next_front;
end

% 保留最好的pareto_front_size个前沿个体
pareto_front_size = size(pareto_front, 1);
if size(new_individuals,1) > pareto_front_size
    [~, ranks] = sortrows(new_individuals(:,num_pareto+1));
    pareto_front_indices = ranks(1:pareto_front_size);
    pareto_front = new_individuals(pareto_front_indices,num_pareto+2:end);
end

% 从帕累托集合中删除被支配的个体
pareto_set_indices = [];
for i = 1:length(pareto_set)
    for j =1:size(pareto_front,1)
        if dominates(calc_fitnesses(pareto_set{i}), pareto_front(j,:))
            pareto_set_indices = [pareto_set_indices, i];
            break
        end
    end
end
pareto_set(pareto_set_indices) = [];

end