

#OLD CODE BELOW HERE -- NEEDS EDITING#

#ready for analysis -- below code rought draft currently

##basic numbers
table(cnee$ratite_accel.1 & cnee$ratite_spec.1, cnee$ratite_accel.2 & cnee$ratite_spec.2)

cnee %>% filter(ratite_accel.2, ratite_spec.2) %>% count(num=floor(ratite_loss_cons_min.prob)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.2, ratite_spec.2) %>% count(num=floor(ratite_loss_cons.prob)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.2, ratite_spec.2) %>% count(num=floor(ratite_loss.prob)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.2, ratite_spec.2) %>% count(num=floor(ratite_loss_cons_min.mat)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.2, ratite_spec.2) %>% count(num=floor(ratite_loss_cons.mat)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.2, ratite_spec.2) %>% count(num=floor(ratite_loss.mat)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))

cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% count(num=floor(ratite_loss_cons_min.prob)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% count(num=floor(ratite_loss_cons.prob)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% count(num=floor(ratite_loss.prob)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% count(num=floor(ratite_loss_cons_min.mat)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% count(num=floor(ratite_loss_cons.mat)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))
cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% count(num=floor(ratite_loss.mat)) %>% summarise(total = sum(n), count = sum(n[num >= 2]), prop = sum(n[num >= 2])/sum(n))

##PLOTS FOR TUFTS TALK##
library(ggthemes)
cnee %>% filter(ratite_accel.1, ratite_spec.1) %>% ggplot(aes(x=floor(ratite_loss_cons.prob))) + geom_bar(fill="red") + coord_flip() + scale_x_reverse() + theme_classic() + theme(text=element_text(size=18))
#raw counts
cnee %>% filter(ratite_accel.2, ratite_spec.2) %>% count(floor(ratite_loss_cons_min.prob))


#neognath distribution of convergent losses
neo_losses <- postaccmat %>%
  mutate(taeGut.loss = (taeGut - taeGut.ficAlb)) %>%
  mutate(ficAlb.loss = (ficAlb - taeGut.pseHum)) %>%
  mutate(pseHum.loss = (pseHum - taeGut.pseHum)) %>%
  mutate(corBra.loss = (corBra - taeGut.corBra)) %>%
  mutate(melUnd.loss = (melUnd - taeGut.melUnd)) %>%
  mutate(falPer.loss = (falPer - taeGut.falPer)) %>%
  mutate(picPub.loss = (picPub - picPub.lepDis)) %>%
  mutate(lepDis.loss = (lepDis - picPub.lepDis)) %>%
  mutate(halLeu.loss = (halLeu - picPub.halLeu)) %>%
  mutate(aptFor.loss = (aptFor - aptFor.pygAde)) %>%
  mutate(pygAde.loss = (pygAde - aptFor.pygAde)) %>%
  mutate(fulGla.loss = (fulGla - aptFor.fulGla)) %>%
  mutate(nipNip.loss = (nipNip - aptFor.nipNip)) %>%
  mutate(balReg.loss = (balReg - balReg.chaVoc)) %>%
  mutate(chaVoc.loss = (chaVoc - balReg.chaVoc)) %>%
  mutate(calAnn.loss = (calAnn - calAnn.chaPel)) %>%
  mutate(chaPel.loss = (chaPel - calAnn.chaPel)) %>%
  mutate(cucCan.loss = (cucCan - calAnn.cucCan)) %>%
  mutate(colLiv.loss = (colLiv - colLiv.mesUni)) %>%
  mutate(mesUni.loss = (mesUni - colLiv.mesUni)) %>% select(taeGut.loss:mesUni.loss)

#sample 10000 tips in batches of 5 and get convergence counts
neo_rand_conv<-data.frame(rep=seq(1,10000), conv_ct_1=NA, conv_ct_2=NA, conv_ct_3=NA, conv_ct_4=NA, conv_ct_5=NA)
for (i in 1:10001) {
  targets <- sample(colnames(neo_losses), 5)
  conv <- neo_losses %>% select(targets) %>% mutate(conv = rowSums(.)) %>% select(conv)
  alllosses <- cbind(neo_losses %>% mutate(all_loss = rowSums(.)) %>% select(all_loss), conv)
  neo_rand_conv$conv_ct_1[i] <- sum(alllosses$conv == 1 & alllosses$all_loss == alllosses$conv) 
  neo_rand_conv$conv_ct_2[i] <- sum(alllosses$conv == 2 & alllosses$all_loss == alllosses$conv)
  neo_rand_conv$conv_ct_3[i] <- sum(alllosses$conv == 3 & alllosses$all_loss == alllosses$conv)
  neo_rand_conv$conv_ct_4[i] <- sum(alllosses$conv == 4 & alllosses$all_loss == alllosses$conv)
  neo_rand_conv$conv_ct_5[i] <- sum(alllosses$conv == 5 & alllosses$all_loss == alllosses$conv)
}

#real ratite data
sum(cnee$ratite_loss_cons.mat == 1 & cnee$tin_loss.mat == 0 & cnee$neo_loss.mat == 0)
sum(cnee$ratite_loss_cons.mat == 2 & cnee$tin_loss.mat == 0 & cnee$neo_loss.mat == 0)
sum(cnee$ratite_loss_cons.mat == 3 & cnee$tin_loss.mat == 0 & cnee$neo_loss.mat == 0)
sum(cnee$ratite_loss_cons.mat == 4 & cnee$tin_loss.mat == 0 & cnee$neo_loss.mat == 0)
sum(cnee$ratite_loss_cons.mat == 5 & cnee$tin_loss.mat == 0 & cnee$neo_loss.mat == 0)

#edit neo rand conv
neo_rand_conv <- neo_rand_conv %>% mutate(conv_tot = conv_ct_2 + conv_ct_3 + conv_ct_4 + conv_ct_5, total_accel = conv_tot + conv_ct_1)
neo_rand_conv <- neo_rand_conv %>% mutate(conv_prop = conv_tot / total_accel)

neo_rand_conv %>% ggplot(aes(x=conv_ct_2)) + geom_histogram(binwidth = 150)
neo_rand_conv %>% ggplot(aes(x=conv_ct_3)) + geom_histogram(binwidth = 10)
neo_rand_conv %>% ggplot(aes(x=conv_ct_4)) + geom_histogram(binwidth = 2)
neo_rand_conv %>% ggplot(aes(x=conv_ct_5)) + geom_histogram(binwidth = 1)

neo_rand_conv %>% ggplot(aes(x=conv_prop)) + geom_histogram(binwidth=0.01, fill="steelblue") + geom_vline(xintercept = (6+118+777+3578)/((6+118+777+3578)+18218), col="red") + theme_classic() + theme(text=element_text(size=18))
neo_rand_conv %>% 
  mutate(conv_prop2 = (conv_ct_3 + conv_ct_4 + conv_ct_5) / (conv_ct_1 + conv_ct_3 + conv_ct_4 + conv_ct_5)) %>% 
  ggplot(aes(x=conv_prop2)) + geom_histogram(binwidth=0.01, fill="steelblue") + geom_vline(xintercept = (6+118+777)/((6+118+777)+18218), col="red") + theme_classic() + theme(text=element_text(size=18))

