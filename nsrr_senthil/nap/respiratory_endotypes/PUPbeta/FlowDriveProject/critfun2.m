function dev = critfun2(X,Y)
model = fitlm(X,Y);
dev = model.Deviance;
end