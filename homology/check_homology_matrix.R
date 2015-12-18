orig<-read.table("original_hog_matrix.txt", header=T)
final<-read.table("updated_hog_matrix.txt", header=T)

orig$total = apply(orig[,2:43], 1, sum)
final$total = apply(final[,2:43], 1, sum)

orig$min = apply(orig[,2:43], 1, min)
final$min = apply(final[,2:43], 1, min)

orig$max = apply(orig[,2:43], 1, max)
final$max = apply(final[,2:43], 1, max)
