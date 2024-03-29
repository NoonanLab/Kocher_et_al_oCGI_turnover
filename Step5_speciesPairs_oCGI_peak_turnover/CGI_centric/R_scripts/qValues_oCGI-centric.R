# 10/17/22
# Purpose: integrate results of permutation test as implemented by speciesPairs_permutation_HPC_CGIcentric.R

setwd('/Users/acadiak/Desktop/CGI/speciesPairs/CGIcentric/permutation')

p_3x3 <- read.table('permutationResults_3x3.txt', header = F)
p_4x3 <- read.table('permutationResults_4x3.txt', header = F)

# add column names
FD_colnames_3x3 <- c('cgiA_peakA.FD','cgiA_peakB.FD','cgiA_peakBoth.FD','cgiB_peakA.FD','cgiB_peakB.FD','cgiB_peakBoth.FD','cgiBoth_peakA.FD','cgiBoth_peakB.FD','cgiBoth_peakBoth.FD')
p_colnames_3x3 <- c('cgiA_peakA.p','cgiA_peakB.p','cgiA_peakBoth.p','cgiB_peakA.p','cgiB_peakB.p','cgiB_peakBoth.p','cgiBoth_peakA.p','cgiBoth_peakB.p','cgiBoth_peakBoth.p')
q_colnames_3x3 <- c('cgiA_peakA.q','cgiA_peakB.q','cgiA_peakBoth.q','cgiB_peakA.q','cgiB_peakB.q','cgiB_peakBoth.q','cgiBoth_peakA.q','cgiBoth_peakB.q','cgiBoth_peakBoth.q')
FD_colnames_4x3 <- c('cgiA_peakA.FD','cgiA_peakB.FD','cgiA_peakBoth.FD','cgiA_peakNeither.FD','cgiB_peakA.FD','cgiB_peakB.FD','cgiB_peakBoth.FD','cgiB_peakNeither.FD','cgiBoth_peakA.FD','cgiBoth_peakB.FD','cgiBoth_peakBoth.FD','cgiBoth_peakNeither.FD')
p_colnames_4x3 <- c('cgiA_peakA.p','cgiA_peakB.p','cgiA_peakBoth.p','cgiA_peakNeither.p','cgiB_peakA.p','cgiB_peakB.p','cgiB_peakBoth.p','cgiB_peakNeither.p','cgiBoth_peakA.p','cgiBoth_peakB.p','cgiBoth_peakBoth.p','cgiBoth_peakNeither.p')
q_colnames_4x3 <- c('cgiA_peakA.q','cgiA_peakB.q','cgiA_peakBoth.q','cgiA_peakNeither.q','cgiB_peakA.q','cgiB_peakB.q','cgiB_peakBoth.q','cgiB_peakNeither.q','cgiBoth_peakA.q','cgiBoth_peakB.q','cgiBoth_peakBoth.q','cgiBoth_peakNeither.q')

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

# current version was run with P = 20000 (see 221017_runPermutationOnHPC.sh)
write.table(q_table_3x3, file='permutation_qAndFD_3x3.txt', quote=F, sep='\t')
write.table(q_table_4x3, file='permutation_qAndFD_4x3.txt', quote=F, sep='\t')
