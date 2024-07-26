function [] = svm_3d_matlab_vis_mod(svmStruct,Xdata,group,justplotsurfaceoncurrentfig)

%originally taken from http://stackoverflow.com/questions/16146212/how-to-plot-a-hyper-plane-in-3d-for-the-svm-results/19969412#19969412
%modified by S. Sands Nov 2016
% modifed a bit more by dwayne mann, 20170630.
% - can handle both the old (svmtrain) and new (fitcsvm) models

switch class(svmStruct)
    case 'ClassificationSVM' % new model
        oldfit = 0;
        sh = svmStruct.Mu; % shift vector
        scalef = svmStruct.Sigma; % scale vector
    case 'struct' % old model
        oldfit = 1;
        sv =  svmStruct.SupportVectors;
        alphaHat = svmStruct.Alpha;
        bias = svmStruct.Bias;
        kfun = svmStruct.KernelFunction;
        kfunargs = svmStruct.KernelFunctionArgs;
        sh = svmStruct.ScaleData.shift; % shift vector
        scalef = svmStruct.ScaleData.scaleFactor; % scale vector
    case 'RegressionLinear' % from fitrlinear
        % not working
        oldfit = 0;
        sh = svmStruct.Bias; % shift vector
        scalef = svmStruct.Beta; % scale vector
    otherwise
        disp('Error in determining model verison');
        return;
end
% set up the data
group = group(~any(isnan(Xdata),2));
Xdata =Xdata(~any(isnan(Xdata),2),:); % remove rows with NaN

% scale and shift data
Xdata1 = repmat(scalef,size(Xdata,1),1).*(Xdata+repmat(sh,size(Xdata,1),1));
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

if oldfit % evaluate the old fit at this point
    f = (feval(kfun,sv,[x y z],kfunargs{:})'*alphaHat(:)) + bias;
end

if iscell(group)
    t = strcmp(group, group{1});
else
    t = group==group(1);
end

% unscale and unshift data
Xdata1 =(Xdata1./repmat(scalef,size(Xdata,1),1)) - repmat(sh,size(Xdata,1),1);
x =(x./repmat(scalef(1),size(x,1),1)) - repmat(sh(1),size(x,1),1);
y =(y./repmat(scalef(2),size(y,1),1)) - repmat(sh(2),size(y,1),1);
z =(z./repmat(scalef(3),size(z,1),1)) - repmat(sh(3),size(z,1),1);

if ~oldfit % evaluate the new fit at this piont
    [~,x_init,~] = predict(svmStruct,[x, y, z]);
    f = x_init(:,1); % if pos class is predictor
    %f = x_init(:,2); % if neg class is predictor
end

if ~justplotsurfaceoncurrentfig
    if oldfit % set them as different figure numbers to allow comparison
        figure(11); clf(figure(11));
    else
        figure(12); clf(figure(12)); 
    end
    plot3(Xdata1(t, 1), Xdata1(t, 2), Xdata1(t, 3), 'b.'); hold on;
    plot3(Xdata1(~t, 1), Xdata1(~t, 2), Xdata1(~t, 3), 'r.');
    if oldfit % only plotting support vectors for old model
        % load unscaled support vectors for plotting
        sv = svmStruct.SupportVectorIndices;
        sv = [Xdata1(sv, :)];
        plot3(sv(:, 1), sv(:, 2), sv(:, 3), 'go');
        title('Old model');
    else
        title('New model');
    end
    
    if iscell(group)
        if oldfit
            legend(group{1},group{end},'support vectors');
        else
            legend(group{1},group{end});
        end
    else
        if oldfit
            legend('group Y','group N','support vectors');
        else
            legend('group Y','group N');
        end
    end
end

x0 = reshape(x, mm);
y0 = reshape(y, mm);
z0 = reshape(z, mm);
v0 = reshape(f, mm);

cutoff = 0;
[faces,verts,~] = isosurface(x0, y0, z0, v0, cutoff, x0); % 0 - continuous variable / cutoff
patch('Vertices', verts, 'Faces', faces, 'FaceColor','k','edgecolor', 'none', 'FaceAlpha', 0.5);

if ~justplotsurfaceoncurrentfig
    grid on
    box on
    legend off
    view(3)
    hold off
end

end