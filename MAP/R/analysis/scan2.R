



load(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
load(file.path(mpath,paste0(pop,'_scan2_bin_mr.rsave')))

plot_test('lod',width=2000,heigh=2000)
plot(bin.em.2,col.scheme = "redblue",contours=T, zlim = c(5,5))
dev.off()

plot_test('lod',width=2000,heigh=2000)
plot(bin.em.2,col.scheme = "redblue", contours=c(16,5), zlim = c(20,18))
dev.off()

plot_test('lod',width=2000,heigh=2000)
plot(bin.imp.2,col.scheme = "redblue", contours=c(16,5), zlim = c(20,18))
dev.off()
