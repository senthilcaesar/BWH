function [] = Plot3Dclassifier(predmodel,Xdata,shift,group,type)

%originally taken from http://stackoverflow.com/questions/16146212/how-to-plot-a-hyper-plane-in-3d-for-the-svm-results/19969412#19969412
%modified by S. Sands Nov 2016

Iincl = sum(isnan(Xdata),2)==0&~isnan(group);

group = group(Iincl);
Xdata = Xdata(Iincl,:); % remove rows with NaN 
% 
% type="svm"; %overwrite

if type=="svm"
    sh = predmodel.ScaleData.shift; % shift vector
    scalef = predmodel.ScaleData.scaleFactor; % scale vector
elseif type=="Quadratic"
    scalef = 1./nanstd(Xdata);
    if isempty(shift)
        sh = -nanmean(Xdata,1);
    else
        sh = shift;%-nanmean(Xdata,1);
    end
end

% scale and shift data
%Xdata1 = repmat(scalef,size(Xdata,1),1).*(Xdata+repmat(sh,size(Xdata,1),1));
Xdata1 = scalef.*(Xdata+sh);

k = 50; 
cubeXMin = min(Xdata1(:,1));
cubeYMin = min(Xdata1(:,2));
cubeZMin = min(Xdata1(:,3));

cubeXMax = max(Xdata1(:,1));
cubeYMax = max(Xdata1(:,2));
cubeZMax = max(Xdata1(:,3));
stepx = (cubeXMax-cubeXMin)/(k-1);
stepy = (cubeYMax-cubeYMin)/(k-1);
stepz = (cubeZMax-cubeZMin)/(k-1);
[x, y, z] = meshgrid(cubeXMin:stepx:cubeXMax,cubeYMin:stepy:cubeYMax,cubeZMin:stepz:cubeZMax);
mm = size(x);
x = x(:);
y = y(:);
z = z(:);

if type=="svm"
sv =  predmodel.SupportVectors;
alphaHat = predmodel.Alpha;
bias = predmodel.Bias;
kfun = predmodel.KernelFunction;
kfunargs = predmodel.KernelFunctionArgs;
f = (feval(kfun,sv,[x y z],kfunargs{:})'*alphaHat(:)) + bias;

% fdata = (feval(kfun,sv,Xdata1,kfunargs{:})'*alphaHat(:)) + bias;
% figure(5)
% plot(fdata,group,'o')

elseif type=="Quadratic"
%myfun = predmodel.model;
myfun = @(A,K,L,Q) K + [A]*L + sum(([A]*Q') .* [A], 2);
K = predmodel.K;
L = predmodel.L;
Q = predmodel.Q;
f=myfun([x y z],K,L,Q);
% fdata=myfun(Xdata1,K,L,Q);
% figure(5)
% plot(fdata,group,'o')
end

% unscale and unshift data 

x = x/scalef(1)-sh(1);
y = y/scalef(2)-sh(2);
z = z/scalef(3)-sh(3);
% 
% 
% 
% Xdata1 =(Xdata1./repmat(scalef,size(Xdata,1),1)) - repmat(sh,size(Xdata,1),1);
if 0
x =(x./repmat(scalef(1),size(x,1),1)) - repmat(sh(1),size(x,1),1);
y =(y./repmat(scalef(2),size(y,1),1)) - repmat(sh(2),size(y,1),1);
z =(z./repmat(scalef(3),size(z,1),1)) - repmat(sh(3),size(z,1),1);
end
x0 = reshape(x, mm);
y0 = reshape(y, mm);
z0 = reshape(z, mm);
v0 = reshape(f, mm);

[faces,verts,colors] = isosurface(x0, y0, z0, v0, 0, x0);
patch('Vertices', verts, 'Faces', faces, 'FaceColor','k','edgecolor', 'none', 'FaceAlpha', 0.5);

