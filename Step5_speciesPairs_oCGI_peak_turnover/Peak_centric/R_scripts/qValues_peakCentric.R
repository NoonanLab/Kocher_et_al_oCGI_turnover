# 1/9/23
# Purpose: integrate results of permutation test as implemented by speciesPairs_permutation_HPC_CGIcentric.R

setwd('/Users/acadiak/Desktop/CGI/speciesPairs/peakCentric/permutation')

p_3x3 <- read.table('permutationResults_3x3.txt', header = F)
p_4x3 <- read.table('permutationResults_4x3.txt', header = F)

# add column names
FD_colnames_3x3 <- c('peakA_cgiA.FD','peakA_cgiB.FD','peakA_cgiBoth.FD','peakB_cgiA.FD','peakB_cgiB.FD','peakB_cgiBoth.FD','peakBoth_cgiA.FD','peakBoth_cgiB.FD','peakBoth_cgiBoth.FD')
p_colnames_3x3 <- c('peakA_cgiA.p','peakA_cgiB.p','peakA_cgiBoth.p','peakB_cgiA.p','peakB_cgiB.p','peakB_cgiBoth.p','peakBoth_cgiA.p','peakBoth_cgiB.p','peakBoth_cgiBoth.p')
q_colnames_3x3 <- c('peakA_cgiA.q','peakA_cgiB.q','peakA_cgiBoth.q','peakB_cgiA.q','peakB_cgiB.q','peakB_cgiBoth.q','peakBoth_cgiA.q','peakBoth_cgiB.q','peakBoth_cgiBoth.q')
FD_colnames_4x3 <- c('peakA_cgiA.FD','peakA_cgiB.FD','peakA_cgiBoth.FD','peakA_cgiNeither.FD','peakB_cgiA.FD','peakB_cgiB.FD','peakB_cgiBoth.FD','peakB_cgiNeither.FD','peakBoth_cgiA.FD','peakBoth_cgiB.FD','peakBoth_cgiBoth.FD','peakBoth_cgiNeither.FD')
p_colnames_4x3 <- c('peakA_cgiA.p','peakA_cgiB.p','peakA_cgiBoth.p','peakA_cgiNeither.p','peakB_cgiA.p','peakB_cgiB.p','peakB_cgiBoth.p','peakB_cgiNeither.p','peakBoth_cgiA.p','peakBoth_cgiB.p','peakBoth_cgiBoth.p','peakBoth_cgiNeither.p')
q_colnames_4x3 <- c('peakA_cgiA.q','peakA_cgiB.q','peakA_cgiBoth.q','peakA_cgiNeither.q','peakB_cgiA.q','peakB_cgiB.q','peakB_cgiBoth.q','peakB_cgiNeither.q','peakBoth_cgiA.q','peakBoth_cgiB.q','peakBoth_cgiBoth.q','peakBoth_cgiNeither.q')

colnames(p_3x3) <- c('speciesPair','speciesA','speciesB','tissue','mark','timePoint', FD_colnames_3x3, p_colnames_3x3)
colnames(p_4x3) <- c('speciesPair','speciesA','speciesB','tissue','mark','timePoint', FD_colnames_4x3, p_colnames_4x3)

# correct p values -> q values
p_matrix_3x3 <- as.matrix(p_3x3[,16:24])
p_adj_vector_3x3 <- p.adjust(p_matrix_3x3, method='BH')
p_adj_3x3 <- matrix(p_adj_vector_3x3, ncol=9)
q_table_3x3 <- p_3x3
q_table_3x3[,16:24] <- p_adj_3x3
q_table_3x3 <- data.frame(q_table_3x3) # contains both FD and q
colnames(q_table_3x3)[16:24] <- q_colnames_3x3

p_matrix_4x3 <- as.matrix(p_4x3[,19:30])
p_adj_vector_4x3 <- p.adjust(p_matrix_4x3,method='BH')
p_adj_4x3 <- matrix(p_adj_vector_4x3, ncol=12)
q_table_4x3 <- p_4x3
q_table_4x3[,19:30] <- p_adj_4x3
q_table_4x3 <- data.frame(q_table_4x3) # contains both FD and q
colnames(q_table_4x3)[19:30] <- q_colnames_4x3

# current version was run with P = 20000 (see 230109_runPermutationOnHPC_peakCentric.sh)
write.table(q_table_3x3, file='permutation_peakCentric_qAndFD_3x3.txt', quote=F, sep='\t')
write.table(q_table_4x3, file='permutation_peakCentric_qAndFD_4x3.txt', quote=F, sep='\t')

