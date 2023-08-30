set <- readRDS("~/Rprojects/ayano/flora_set.rds")


AS_d1 = set$points %>% filter(group == "PL",ENA_DIRECTION == "response") %>% as.matrix()
AS_d1 = AS_d1[,1]

AS_d2 = set$points %>% filter(group == "PL",ENA_DIRECTION == "response") %>% as.matrix()
AS_d2 = AS_d2[,2]


CN_d1 = set$points %>% filter(group == "CN",ENA_DIRECTION == "response") %>% as.matrix()
CN_d1 = CN_d1[,1]

CN_d2 = set$points %>% filter(group == "CN",ENA_DIRECTION == "response") %>% as.matrix()
CN_d2 = CN_d2[,2]

FS_d2 = set$points %>% filter(group == "GE",ENA_DIRECTION == "response") %>% as.matrix()
FS_d2 = FS_d2[,2]



#####################


AS_d1 <- as.matrix(set$points$group$PL)[,1]
AS_d2 <- as.matrix(set$points$group$PL)[,2]

CN_d1 <- as.matrix(set$points$group$CN)[,1]
CN_d2 <- as.matrix(set$points$group$CN)[,2]

FS_d1 <- as.matrix(set$points$group$GE)[,1]
FS_d2 <- as.matrix(set$points$group$GE)[,2]

AS_S4_d1 <- as.matrix(set$points[group == "PL" & scaffolding == "S4"])[,1]
AS_S4_d2 <- as.matrix(set$points[group == "PL" & scaffolding == "S4"])[,2]

AS_S5_d1 <- as.matrix(set$points[group == "PL" & scaffolding == "S5"])[,1]
AS_S5_d2 <- as.matrix(set$points[group == "PL" & scaffolding == "S5"])[,2]

CN_S4_d1 <- as.matrix(set$points[group == "CN" & scaffolding == "S4"])[,1]
CN_S4_d2 <- as.matrix(set$points[group == "CN" & scaffolding == "S4"])[,2]

CN_S5_d1 <- as.matrix(set$points[group == "CN" & scaffolding == "S5"])[,1]
CN_S5_d2 <- as.matrix(set$points[group == "CN" & scaffolding == "S5"])[,2]

FS_S4_d1 <- as.matrix(set$points[group == "GE" & scaffolding == "S4"])[,1]
FS_S4_d2 <- as.matrix(set$points[group == "GE" & scaffolding == "S4"])[,2]

FS_S5_d1 <- as.matrix(set$points[group == "GE" & scaffolding == "S5"])[,1]
FS_S5_d2 <- as.matrix(set$points[group == "GE" & scaffolding == "S5"])[,2]

# 5. Tukey's Test 

## 5.1 task level 

AS_d1 <- as.data.frame(AS_d1)
AS_d1$group <- "AS_d1"
names(AS_d1)[1] <- "pts"

CN_d1 <- as.data.frame(CN_d1)
CN_d1$group <- "CN_d1"
names(CN_d1)[1] <- "pts"

FS_d1 <- as.data.frame(FS_d1)
FS_d1$group <- "FS_d1"
names(FS_d1)[1] <- "pts"

task <- data.table::rbindlist(list(AS_d1, CN_d1, FS_d1))

lm <- lm(pts ~ group, data = task)
aov <- aov(lm)
summary(aov)

tukey.test <- TukeyHSD(aov)
tukey.test

plot(tukey.test)

AS_d2 <- as.data.frame(AS_d2)
AS_d2$group <- "AS_d2"
names(AS_d2)[1] <- "pts"

CN_d2 <- as.data.frame(CN_d2)
CN_d2$group <- "CN_d2"
names(CN_d2)[1] <- "pts"

FS_d2 <- as.data.frame(FS_d2)
FS_d2$group <- "FS_d2"
names(FS_d2)[1] <- "pts"

task <- data.table::rbindlist(list(AS_d2, CN_d2, FS_d2))

lm <- lm(pts ~ group, data = task)
aov <- aov(lm)
summary(aov)

tukey.test <- TukeyHSD(aov)
tukey.test

plot(tukey.test)


