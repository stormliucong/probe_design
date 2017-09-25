library(reshape2)
mat <- read.table("SEQ4416.flag.stat",head=F)
mat2 <-dcast(mat,V1~V2)

# 69 = read paired
# read unmapped
# first in pair

# 73=read paired
# mate unmapped
# first in pair

# 77 = read paired
# read unmapped
# mate unmapped
# first in pair



pair_map_r1 <- mat2$`99`
pair_map_r1_reverse <- mat2$`83`
pair_map_r2 <- mat2$`147`
pair_map_r2_reverse <- mat2$`163`


r2_discorcondant <- mat2$`97`
r2_discorcondant2 <- mat2$`113`
r1_discorcondant <- mat2$`145`
r1_discorcondant_2 <- mat2$`161`
r1_discorcondant_3 <- mat2$`177`



names(mate_discorcondant) <- mat2$V1


matrix <- diag(pair_map_r1) + diag(pair_map_r2)
colnames(matrix) <- mat2$V1
rownames(matrix) <- mat2$V1


