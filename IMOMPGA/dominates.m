function result = dominates(fitness1, fitness2)
% 判断一个个体是否支配另一个个体

result = all(fitness1 <= fitness2) && any(fitness1 < fitness2);

end