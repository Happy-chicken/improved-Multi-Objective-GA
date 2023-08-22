# 用改进的多种群遗传算法解决车间调度问题
# Paper title: Fuzzy Target Job Shop Scheduling Based on Improved Multi-population Genetic Algorithm with Clustering
abstract: Owing to the traditional genetic algorithm is difficult to hit the optimal solution in the optimization problem of fuzzy target shop scheduling, this paper converts the weighting of the two goals of processing completion time and customer satisfaction into a single objective task, and adds fuzzy completion time constraints on the basis of traditional shop scheduling constraints A fuzzy objective shop scheduling model is constructed, some genetic operators using different strategies are designed and a multi-population genetic algorithm based on clustering is proposed according to the similarity of populations The algorithm first generates multiple initial populations by extending the integer encoding of the processing order and processing machines, and each population undergoes independent genetic operations and then forms importation, integration, and artistic selection to share the excellent genes of individuals The algorithm is verified through experimental simulations that it has strong stability and global search capability: The probability of finding the global optimal solution has increased to 61%.
# file description
GA.m: 普通遗传算法主函数
MPGA.m: 多种群遗传算法主函数
IMPGA.m: 利用K-means聚类改进后的多种群遗传算法主函数
