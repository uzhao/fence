load("X.add.n.75.rda")
load("X.new.1.n.75.rda")

beta = c(1,1,1,1,.1,.05,.25)

n = nrow(X.new.1)
p = ncol(X.new.1)


M.1 = 1:7
M.2 = c(8:15)
M.3 = c(16:22)
M.4 = c(23:30)

m.ind = 1:30

template = function(M.ind,X,z,z.z,z.b,z.z.b){

M.ind.2 = m.ind[-M.ind]	
X.1 = X[,M.ind]
X.2 = X[,M.ind.2]
X.1.1 = t(X.1)%*%X.1
X.1.2 = t(X.1)%*%X.2
X.2.2 = t(X.2)%*%X.2
X.2.2.inv = solve(X.2.2)
X.1.z = t(X.1)%*%z
X.2.z = t(X.2)%*%z
model.cand = model.dic(M.ind)
Q.est = sapply(model.cand,Q.hat,X,X.2,z,z.z,X.2.z,X.2.2,X.2.2.inv)
len.Q = length(Q.est)
Q.diff = Q.est - Q.est[len.Q]
X.2.z.b = apply(z.b,2,X.2.z.boot,X.2=X.2)

Q.diff.b = matrix(0,nrow = len.Q, ncol = 100)

for(k in 1:100){
	Q.est.b = sapply(model.cand,Q.hat,X,X.2,z.b[,k],z.z.b[k],X.2.z.b[,k],X.2.2,X.2.2.inv)		
	Q.diff.b[,k] = Q.est.b - Q.est.b[len.Q]
}

upper.Cn.1 = max(Q.diff.b)
lower.Cn.1 = min(Q.diff.b)
Cn = seq(lower.Cn.1,upper.Cn.1,by=10)
len.Cn = length(Cn)
p.star = numeric(len.Cn)
M.index = matrix(0,nrow=100,ncol=len.Cn)


for(k in 1:100){
	for(w in 1:len.Cn){
		M.index[k,w] = model.sel(model.cand,Q.diff.b[,k],Cn[w])
	}
}

for(j in 1:len.Cn){
	p.star[j] =prob.ind(M.index[,j])
}


p.star.adj = p.star[-1]
Cn.adj = Cn[-1]
len.Cn.adj = length(Cn.adj)

#test = loess(p.star.adj~Cn.adj,span=.02,degree=1)
#temp.s = round(test$fitted, digits = 2)
#lines(Cn.adj,temp.s,pch=20,lwd=2,col="red")
temp.s = p.star.adj


Cn.peak = find.peak(temp.s,len.Cn.adj)
len.Cn.peak = length(Cn.peak)

if(len.Cn.peak == 0){
	max.p.index = p.index(temp.s)
	Cn.range = Cn.adj[max.p.index]
	Cn.star= median(Cn.range)
}else if(len.Cn.peak == 1){
	a.star = find.a.star(temp.s,len.Cn.adj)
	b.star = find.b.star(temp.s,len.Cn.adj,a.star)
	range.n = a.star:b.star
	temp.1 = temp.s[range.n]
	Cn.1 = Cn.adj[range.n]
	max.p.index = p.index(temp.1)
	Cn.range = Cn.1[max.p.index]
	Cn.star = median(Cn.range)
}else{
	Cn.star = Cn.adj[Cn.peak[1]]
	p.star.peak = temp.s[Cn.peak]
	p.hat = max(p.star.peak)
	#Cn.star = Cn.adj[Cn.peak[which(p.star.peak == p.hat)][1]]
}
index.m = model.sel(model.cand,Q.diff,Cn.star)
if(index.m == 0){
	index.m = model.c(model.cand,Q.diff,len.Q,Cn.star)
}
list(model.cand[[index.m]],p.star)
}

template.2 = function(M.ind,X,z,z.z,z.b,z.z.b){

M.ind.2 = m.ind[-M.ind]	
X.1 = X[,M.ind]
X.2 = X[,M.ind.2]
X.1.1 = t(X.1)%*%X.1
X.1.2 = t(X.1)%*%X.2
X.2.2 = t(X.2)%*%X.2
X.2.2.inv = solve(X.2.2)
X.1.z = t(X.1)%*%z
X.2.z = t(X.2)%*%z
model.cand = model.dic(M.ind)
Q.est = sapply(model.cand,Q.hat,X,X.2,z,z.z,X.2.z,X.2.2,X.2.2.inv)
len.Q = length(Q.est)
Q.diff = Q.est - Q.est[len.Q]
X.2.z.b = apply(z.b,2,X.2.z.boot,X.2=X.2)

Q.diff.b = matrix(0,nrow = len.Q, ncol = 100)

for(k in 1:100){
	Q.est.b = sapply(model.cand,Q.hat,X,X.2,z.b[,k],z.z.b[k],X.2.z.b[,k],X.2.2,X.2.2.inv)		
	Q.diff.b[,k] = Q.est.b - Q.est.b[len.Q]
}

upper.Cn.1 = max(Q.diff.b)
lower.Cn.1 = min(Q.diff.b)
Cn = seq(lower.Cn.1,upper.Cn.1,by=1)
len.Cn = length(Cn)
p.star = numeric(len.Cn)
M.index = matrix(0,nrow=100,ncol=len.Cn)


for(k in 1:100){
	for(w in 1:len.Cn){
		M.index[k,w] = model.sel(model.cand,Q.diff.b[,k],Cn[w])
	}
}

for(j in 1:len.Cn){
	p.star[j] =prob.ind(M.index[,j])
}


p.star.adj = p.star[-c(1:15)]
Cn.adj = Cn[-(1:15)]
len.Cn.adj = length(Cn.adj)

test = loess(p.star.adj~Cn.adj,span=.1,degree=1)
temp.s = round(test$fitted, digits = 2)
#lines(Cn.adj,temp.s,pch=20,lwd=2,col="red")


Cn.peak = find.peak(temp.s,len.Cn.adj)
len.Cn.peak = length(Cn.peak)

if(len.Cn.peak == 0){
	max.p.index = p.index(temp.s)
	Cn.range = Cn.adj[max.p.index]
	Cn.star= median(Cn.range)
}else if(len.Cn.peak == 1){
	a.star = find.a.star(temp.s,len.Cn.adj)
	b.star = find.b.star(temp.s,len.Cn.adj,a.star)
	range.n = a.star:b.star
	temp.1 = temp.s[range.n]
	Cn.1 = Cn.adj[range.n]
	max.p.index = p.index(temp.1)
	Cn.range = Cn.1[max.p.index]
	Cn.star = median(Cn.range)
}else{
	Cn.star = Cn.adj[Cn.peak[1]]
	p.star.peak = temp.s[Cn.peak]
	p.hat = max(p.star.peak)
	#Cn.star = Cn.adj[Cn.peak[which(p.star.peak == p.hat)][1]]
}
index.m = model.sel(model.cand,Q.diff,Cn.star)
if(index.m == 0){
	index.m = model.c(model.cand,Q.diff,len.Q,Cn.star)
}
list(model.cand[[index.m]],p.star)
}


template.3 = function(M.ind,X,z,z.z,z.b,z.z.b){

M.ind.2 = m.ind[-M.ind]	
X.1 = X[,M.ind]
X.2 = X[,M.ind.2]
X.1.1 = t(X.1)%*%X.1
X.1.2 = t(X.1)%*%X.2
X.2.2 = t(X.2)%*%X.2
X.2.2.inv = solve(X.2.2)
X.1.z = t(X.1)%*%z
X.2.z = t(X.2)%*%z
model.cand = model.dic(M.ind)
Q.est = sapply(model.cand,Q.hat,X,X.2,z,z.z,X.2.z,X.2.2,X.2.2.inv)
len.Q = length(Q.est)
Q.diff = Q.est - Q.est[len.Q]
X.2.z.b = apply(z.b,2,X.2.z.boot,X.2=X.2)

Q.diff.b = matrix(0,nrow = len.Q, ncol = 100)

for(k in 1:100){
	Q.est.b = sapply(model.cand,Q.hat,X,X.2,z.b[,k],z.z.b[k],X.2.z.b[,k],X.2.2,X.2.2.inv)		
	Q.diff.b[,k] = Q.est.b - Q.est.b[len.Q]
}

upper.Cn.1 = max(Q.diff.b)
lower.Cn.1 = min(Q.diff.b)
Cn = seq(lower.Cn.1,upper.Cn.1,by=1)
len.Cn = length(Cn)
p.star = numeric(len.Cn)
M.index = matrix(0,nrow=100,ncol=len.Cn)


for(k in 1:100){
	for(w in 1:len.Cn){
		M.index[k,w] = model.sel(model.cand,Q.diff.b[,k],Cn[w])
	}
}

for(j in 1:len.Cn){
	p.star[j] =prob.ind(M.index[,j])
}


p.star.adj = p.star[-1]
Cn.adj = Cn[-1]
len.Cn.adj = length(Cn.adj)

test = loess(p.star.adj~Cn.adj,span=.5,degree=1)
temp.s = round(test$fitted, digits = 2)
#lines(Cn.adj,temp.s,pch=20,lwd=2,col="red")

Cn.peak = find.peak(temp.s,len.Cn.adj)
len.Cn.peak = length(Cn.peak)

if(len.Cn.peak == 0){
	max.p.index = p.index(temp.s)
	Cn.range = Cn.adj[max.p.index]
	Cn.star= median(Cn.range)
}else if(len.Cn.peak == 1){
	a.star = find.a.star(temp.s,len.Cn.adj)
	b.star = find.b.star(temp.s,len.Cn.adj,a.star)
	range.n = a.star:b.star
	temp.1 = temp.s[range.n]
	Cn.1 = Cn.adj[range.n]
	max.p.index = p.index(temp.1)
	Cn.range = Cn.1[max.p.index]
	Cn.star = median(Cn.range)
}else{
	#Cn.star = Cn.adj[Cn.peak[1]]
	p.star.peak = temp.s[Cn.peak]
	p.hat = max(p.star.peak)
	Cn.star = Cn.adj[Cn.peak[which(p.star.peak == p.hat)][1]]
}
index.m = model.sel(model.cand,Q.diff,Cn.star)
if(index.m == 0){
	index.m = model.c(model.cand,Q.diff,len.Q,Cn.star)
}
list(model.cand[[index.m]],p.star)
}

template.4 = function(M.ind,X,z,z.z,z.b,z.z.b){

M.ind.2 = m.ind[-M.ind]	
X.1 = X[,M.ind]
X.2 = X[,M.ind.2]
X.1.1 = t(X.1)%*%X.1
X.1.2 = t(X.1)%*%X.2
X.2.2 = t(X.2)%*%X.2
X.2.2.inv = solve(X.2.2)
X.1.z = t(X.1)%*%z
X.2.z = t(X.2)%*%z
model.cand = model.dic(M.ind)
Q.est = sapply(model.cand,Q.hat,X,X.2,z,z.z,X.2.z,X.2.2,X.2.2.inv)
len.Q = length(Q.est)
Q.diff = Q.est - Q.est[len.Q]
X.2.z.b = apply(z.b,2,X.2.z.boot,X.2=X.2)

Q.diff.b = matrix(0,nrow = len.Q, ncol = 100)

for(k in 1:100){
	Q.est.b = sapply(model.cand,Q.hat,X,X.2,z.b[,k],z.z.b[k],X.2.z.b[,k],X.2.2,X.2.2.inv)		
	Q.diff.b[,k] = Q.est.b - Q.est.b[len.Q]
}

upper.Cn.1 = max(Q.diff.b)
lower.Cn.1 = min(Q.diff.b)
Cn = seq(0,200,by=1)
len.Cn = length(Cn)
p.star = numeric(len.Cn)
M.index = matrix(0,nrow=100,ncol=len.Cn)


for(k in 1:100){
	for(w in 1:len.Cn){
		M.index[k,w] = model.sel(model.cand,Q.diff.b[,k],Cn[w])
	}
}

for(j in 1:len.Cn){
	p.star[j] =prob.ind(M.index[,j])
}


p.star.adj = p.star[-c(1:15)]
Cn.adj = Cn[-(1:15)]
len.Cn.adj = length(Cn.adj)

test = loess(p.star.adj~Cn.adj,span=.1,degree=1)
temp.s = round(test$fitted, digits = 2)
#lines(Cn.adj,temp.s,pch=20,lwd=2,col="red")


Cn.peak = find.peak(temp.s,len.Cn.adj)
len.Cn.peak = length(Cn.peak)

if(len.Cn.peak == 0){
	max.p.index = p.index(temp.s)
	Cn.range = Cn.adj[max.p.index]
	Cn.star= median(Cn.range)
}else if(len.Cn.peak == 1){
	a.star = find.a.star(temp.s,len.Cn.adj)
	b.star = find.b.star(temp.s,len.Cn.adj,a.star)
	range.n = a.star:b.star
	temp.1 = temp.s[range.n]
	Cn.1 = Cn.adj[range.n]
	max.p.index = p.index(temp.1)
	Cn.range = Cn.1[max.p.index]
	Cn.star = median(Cn.range)
}else{
	Cn.star = Cn.adj[Cn.peak[1]]
	p.star.peak = temp.s[Cn.peak]
	p.hat = max(p.star.peak)
	#Cn.star = Cn.adj[Cn.peak[which(p.star.peak == p.hat)][1]]
}
index.m = model.sel(model.cand,Q.diff,Cn.star)
if(index.m == 0){
	index.m = model.c(model.cand,Q.diff,len.Q,Cn.star)
}
list(model.cand[[index.m]],p.star)
}

model.dic = function(m.sel){
len.m = length(m.sel)
num = numeric(len.m)
index.1 = vector("list", len.m)

for(i in 1:len.m){
index.1[[i]] = combn(m.sel,i)
num[i] = ncol(index.1[[i]])
}

tot.n = 1 + sum(num)
model = vector("list", tot.n)

model[[1]] = NULL
a = 1
for(i in 1:len.m){
	for(j in 1:num[i]){
		a = a + 1
		model[[a]] = index.1[[i]][,j]
	}
}
model
}

Q.hat=function(M,X,X.2,z,z.z,X.2.z,X.2.2,X.2.2.inv){

	if(length(M) == 0){
		QM.part.1 = z.z 
		QM.part.2 = t(X.2.z)%*%X.2.2.inv%*%X.2.z
	}else{
		QM.part.1 = z.z - t(X.2.z)%*%X.2.2.inv%*%X.2.z
		X.1.s = X[,M]

		X.1.1.s = t(X.1.s)%*%X.1.s
		X.1.2.s = t(X.1.s)%*%X.2

		X.1.z.s = t(X.1.s)%*%z

		part.1.b.s = X.1.1.s - X.1.2.s%*%X.2.2.inv%*%t(X.1.2.s)
		part.2.b.s = X.1.z.s - X.1.2.s%*%X.2.2.inv%*%X.2.z
		
		QM.part.2 = t(part.2.b.s)%*%solve(part.1.b.s)%*%part.2.b.s

	}
as.numeric(QM.part.1 - QM.part.2)
}

X.2.z.boot = function(A,X.2){
t(X.2)%*%A
}

dim.m = function(X){
length(X)
}

model.sel = function(model.cand,Q.diff.star,c){
	len.dim = sapply(model.cand,dim.m)
	m.index = 0
	uniq.dim = unique(len.dim)
	len.uniq.dim = length(uniq.dim)
	j = 1
	while((j != 0) && (j <= len.uniq.dim )){
		dim.M = uniq.dim[j]
		index = which(len.dim == dim.M)
		Q.diff.sub = Q.diff.star[index]
		pos = Q.diff.sub <= c
		index.s = which(pos == 1)
		if(length(index.s) == 0){
			j = j +1
		}else{
			m.index = which(Q.diff.star == min(Q.diff.sub[index.s]))[1]
			j = 0
		}
	}
m.index
}

model.c = function(model.cand,Q.diff.star,len.Q,c){
	diff = numeric(len.Q)
	for(j in 1:len.Q){
		diff[j] = Q.diff.star[j] - c
	}
	m.index = which(diff == min(diff))
m.index
}


prob.ind = function(X){
prob = 0
a = table(X)
b = max(a)
ind = as.numeric(names(which(a == b)))[1]
if(ind == 0){
	prob = 0
}else{
	prob = b/length(X)
}
prob
}

find.peak = function(temp,len.Cn){
index = NULL

a = NULL
j = 2
	if(temp[(j-1)] < temp[j]){
		a= j-1
	}else{
		while((j!=0) &&(j <= (len.Cn-1))){
			if((temp[j]<= temp[(j-1)]) && (temp[j] < temp[(j+1)])){
				a = j
				j =0
			}else{
				j = j+1
			}
		}	
	}
if(length(a) == 0){
a = 0
b = NULL
j = a+1
u = 1

while(u != 0){
	
	while((j!=0)&&(j <= (len.Cn-1))){
		if((temp[j]< temp[(j+1)]) && (temp[j] <= temp[(j-1)])){
			b = j 
			j =0
		}else{
			j = j+1
		}
	}
	
	if(length(b) == 0){
		u = 0
	}else{
		range.s = a:b
		temp.s = temp[range.s]
		max.p = max(temp.s)
		index = c(index,floor(median(range.s[which(temp.s == max.p)])))
	
		j = b 

		while((j!=0) &&(j <= (len.Cn-1))){
			if((temp[j]<= temp[(j-1)]) && (temp[j] < temp[(j+1)])){
				a = j
				j =0
			}else{
				j = j+1
				if(j == len.Cn){
					u = 0
				}
			}	
		}
	
		b = NULL
		j = a +1

	}
}	
}else{
b = NULL
j = a+1
u = 1

while(u != 0){
	
	while((j!=0)&&(j <= (len.Cn-1))){
		if((temp[j]<= temp[(j+1)]) && (temp[j] < temp[(j-1)])){
			b = j 
			j =0
		}else{
			j = j+1
		}
	}
	
	if(length(b) == 0){
		u = 0
	}else{
		range.s = a:b
		temp.s = temp[range.s]
		max.p = max(temp.s)
		index = c(index,floor(median(range.s[which(temp.s == max.p)])))
	
		j = b 

		while((j!=0) &&(j <= (len.Cn-1))){
			if((temp[j]<= temp[(j-1)]) && (temp[j] < temp[(j+1)])){
				a = j
				j =0
			}else{
				j = j+1
				if(j == len.Cn){
					u = 0
				}
			}	
		}
	
		b = NULL
		j = a +1

	}	
}
}
index
}

find.a.star = function(temp,len.Cn){
a.star = NULL
j = 2

if(temp[(j-1)] < temp[j]){
	a.star = j-1
}else{
	while((j!=0) &&(j <= (len.Cn-1))){
		if((temp[j]<= temp[(j-1)]) && (temp[j] < temp[(j+1)])){
			a.star = j
			j =0
		}else{
			j = j+1
		}	
	}
}
a.star
}

find.b.star = function(temp,len.Cn,a.star){

b.star = NULL
j =len.Cn 
if(length(a.star)>0){	
	if(temp[j] < temp[(j-1)]){
		b.star = j
		j = 0
	}else{
		j = j - 1
		while((j!=0) && (j > a.star)){
			if((temp[j]<= temp[(j+1)]) && (temp[j] < temp[(j-1)])){
				b.star = j 
				j =0
			}else{
				j = j-1
			}
		}		
	}
}else{ 
	if(temp[j] < temp[(j-1)]){
		b.star = j
		j = 0
	}else{
		j = j-1
		while((j!=0) &&(j > 1)){
			if((temp[j]<= temp[(j+1)]) && (temp[j] < temp[(j-1)])){
				b.star = j 
				j =0
			}else{
				j = j-1
			}	
		}
	}
}
b.star
}

p.index = function(X){
max.p.index = which(X == max(X))
sub.f = max.p.index[-length(max.p.index)]
sub.l = max.p.index[-1]
diff = sub.l - sub.f
index = which(diff != 1)
len.index = length(index)
if(len.index == 0){
max.p.index = max.p.index
}else{
cut.off = index[1] 
max.p.index = max.p.index[1:cut.off]
}
max.p.index
}

Q.add.hat = function(M,X.t,z){
X.sub = X.t[,M]
b.hat = solve(t(X.sub)%*%X.sub)%*%t(X.sub)%*%z
t(z - X.sub%*%b.hat)%*%(z - X.sub%*%b.hat)
}

model.dic.2 = function(m.sel){
len.m = length(m.sel)
num = numeric(len.m)
index.1 = vector("list", len.m)

for(i in 1:len.m){
index.1[[i]] = combn(m.sel,i)
num[i] = ncol(index.1[[i]])
}

tot.n = sum(num)
model = vector("list", tot.n)

a = 0
for(i in 1:len.m){
	for(j in 1:num[i]){
		a = a + 1
		model[[a]] = index.1[[i]][,j]
	}
}
model
}

result = vector("list", 100)
for(r in 1:100){
set.seed(r)
v_i = rep(rnorm(n/3),each = 3)
e  = rnorm(n)
X.t = cbind(X.new.1[,c(1:4)],X.new.1[,25],X.new.1[,c(7:8)])
y = X.t %*% beta + v_i + e
y.y = t(y)%*%y
beta.hat = solve(t(X.new.1)%*%X.new.1)%*%t(X.new.1)%*%y
eps.b = matrix(rnorm(n*100, mean = 0, sd =1),nrow=n)
y.b = matrix(0,nrow= n,ncol=100)
y.y.b = numeric(100)
for(i in 1:100){
y.b[,i] = X.new.1%*%beta.hat + eps.b[,i]
y.y.b[i] = t(y.b[,i])%*%y.b[,i]
}

t.1 = template(M.1,X.new.1,y,y.y,y.b,y.y.b)
mar.sel.1 = t.1[[1]]
pstar.1 = t.1[[2]]
t.2 = template.2(M.2,X.new.1,y,y.y,y.b,y.y.b)
mar.sel.2 = t.2[[1]]
pstar.2 = t.2[[2]]
t.3 = template.3(M.3,X.new.1,y,y.y,y.b,y.y.b)
mar.sel.3 = t.3[[1]]
pstar.3 = t.3[[2]]
t.4 = template.4(M.4,X.new.1,y,y.y,y.b,y.y.b)
mar.sel.4 = t.4[[1]]
pstar.4 = t.4[[2]]
p.star = list(pstar.1,pstar.2,pstar.3,pstar.4)

marker.sel = c(mar.sel.1,mar.sel.2,mar.sel.3,mar.sel.4)

result.1 = list(marker.sel,p.star)
X.base = cbind(X.new.1,X.add)
m.1 = c(marker.sel,ncol(X.base))
Q.f.add = Q.add.hat(m.1,X.base,y)
model.cand.2 = model.dic.2(marker.sel)
len.cand.2 = length(model.cand.2)
Q.2.est = sapply(model.cand.2,Q.add.hat,X.t=X.new.1,z=y)
Q.diff.2 = Q.2.est - Q.f.add
Q.diff.2.b = matrix(0,nrow=len.cand.2,ncol=100)
for(i in 1:100){
	temp.1 = sapply(model.cand.2,Q.add.hat,X.t=X.new.1,z=y.b[,i])
	temp.2 = Q.add.hat(m.1,X.base,y.b[,i])
	Q.diff.2.b[,i] = temp.1 - temp.2 
}

upper.Cn.2 = max(Q.diff.2.b)
Cn.2.1 = seq(0,200,by=1)
Cn.2.2 = seq(201,(upper.Cn.2 - 50),by=10000)
Cn.2.3 = seq((upper.Cn.2 - 49),upper.Cn.2,by=1)
Cn.2 = c(Cn.2.1,Cn.2.2,Cn.2.3)
#Cn.2 = seq(0,upper.Cn.2,by=1000)
len.Cn.2 = length(Cn.2)

p.star.2 = numeric(len.Cn.2)
M.index.2 = matrix(0,nrow=100,ncol=len.Cn.2)

for(k in 1:100){
	for(w in 1:len.Cn.2){
		M.index.2[k,w] = model.sel(model.cand.2,Q.diff.2.b[,k],Cn.2[w])
	}
}

for(j in 1:len.Cn.2){
	p.star.2[j] =prob.ind(M.index.2[,j])
}
p.star.2.adj = p.star.2[-c(1:5)]
Cn.2.adj = Cn.2[-c(1:5)]
len.Cn.2.adj = length(Cn.2.adj)

temp.4 = p.star.2.adj
Cn.peak.2 = find.peak(temp.4,len.Cn.2.adj)
len.Cn.peak.2 = length(Cn.peak.2)

if(len.Cn.peak.2 == 0){
	max.p.index.2 = p.index(temp.4)
	Cn.range.2 = Cn.2.adj[max.p.index.2]
	Cn.star.2 = median(Cn.range.2)
	Cn.star.2.1 = Cn.star.2
	Cn.star.2.2 = Cn.star.2
	Cn.star.2.3 = Cn.star.2
	Cn.star.2.4 = Cn.star.2	
}else if(len.Cn.peak.2 == 1){
	a.star.2 = find.a.star(temp.4,len.Cn.2.adj)
	b.star.2 = find.b.star(temp.4,len.Cn.2.adj,a.star.2)
	range.n.2 = a.star.2:b.star.2
	temp.1.2 = temp.4[range.n.2]
	Cn.1.2 = Cn.2.adj[range.n.2]
	max.p.index.2 = p.index(temp.1.2)
	Cn.range.2 = Cn.1.2[max.p.index.2]
	Cn.star.2 = median(Cn.range.2)
	Cn.star.2.1 = Cn.star.2
	Cn.star.2.2 = Cn.star.2
	Cn.star.2.3 = Cn.star.2
	Cn.star.2.4 = Cn.star.2	
}else{
	p.star.peak.2 = temp.4[Cn.peak.2]
	p.hat.2 = max(p.star.peak.2)
	Cn.star.2 = Cn.2.adj[Cn.peak.2[which(p.star.peak.2 == p.hat.2)][1]]
	if(p.hat.2 == 1){
		p.hat.2 = p.hat.2 - 10^(-4)
	}
	lower = p.hat.2 - 1.96*sqrt((p.hat.2*(1-p.hat.2))/100)
	lower.1 = p.hat.2 - 2.33*sqrt((p.hat.2*(1-p.hat.2))/100)
	lower.2 = p.hat.2 - 2.58*sqrt((p.hat.2*(1-p.hat.2))/100)

	Cn.star.2.1 = Cn.2.adj[Cn.peak.2[which(p.star.peak.2 >= lower)][1]]
	Cn.star.2.2 = Cn.2.adj[Cn.peak.2[which(p.star.peak.2 >= lower.1)][1]]
	Cn.star.2.3 = Cn.2.adj[Cn.peak.2[which(p.star.peak.2 >= lower.2)][1]]
	Cn.star.2.4 = Cn.2.adj[Cn.peak.2[1]]
}

index.m.2 = model.sel(model.cand.2,Q.diff.2,Cn.star.2)
if(index.m.2 == 0){
	index.m.2 = model.c(model.cand.2,Q.diff.2,len.cand.2,Cn.star.2)
}

index.m.2.1 = model.sel(model.cand.2,Q.diff.2,Cn.star.2.1)
if(index.m.2.1 == 0){
	index.m.2.1 = model.c(model.cand.2,Q.diff.2,len.cand.2,Cn.star.2.1)
}

index.m.2.2 = model.sel(model.cand.2,Q.diff.2,Cn.star.2.2)
if(index.m.2.2 == 0){
	index.m.2.2 = model.c(model.cand.2,Q.diff.2,len.cand.2,Cn.star.2.2)
}

index.m.2.3 = model.sel(model.cand.2,Q.diff.2,Cn.star.2.3)
if(index.m.2.3 == 0){
	index.m.2.3 = model.c(model.cand.2,Q.diff.2,len.cand.2,Cn.star.2.3)
}

index.m.2.4 = model.sel(model.cand.2,Q.diff.2,Cn.star.2.4)
if(index.m.2.4 == 0){
	index.m.2.4 = model.c(model.cand.2,Q.diff.2,len.cand.2,Cn.star.2.4)
}

mar.sel.2 = model.cand.2[[index.m.2]]
mar.sel.2.1 = model.cand.2[[index.m.2.1]]
mar.sel.2.2 = model.cand.2[[index.m.2.2]]
mar.sel.2.3 = model.cand.2[[index.m.2.3]]
mar.sel.2.4 = model.cand.2[[index.m.2.4]]

mar.sel.s = list(mar.sel.2,mar.sel.2.1,mar.sel.2.2,mar.sel.2.3,mar.sel.2.4)
Cn.star.s = c(Cn.star.2,Cn.star.2.1,Cn.star.2.2,Cn.star.2.3,Cn.star.2.4)
result.2 = list(mar.sel.s,Cn.star.s,p.star.2)
result[[r]] = list(result.1,result.2)
save(result, file="result.RF.LDA.n.150.rda")
}




