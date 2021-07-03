function files_out = spmup_hrf_boost_yc(varargin)
% 2020.04.20 some functions fixed 
% this function takes a SPM.mat and related files to create boosted beta
% parameters and con images if the hrf doesn't fit properly the data, some
% information can be recovered using the 1st and 2nd derivatives - most of
% the time this boosts the hrf value but it can also reduce it if there are
% strong shape differences - it is thus recommended to use 1st derivative
% only, unless you have a slow design.
%
% Ref.  Pernet CR (2014) Misconceptions in the use of the General Linear
% Model applied to functional MRI: a tutorial for junior neuro-imagers.
% Front. Neurosci. 8:1. doi: 10.3389/fnins.2014.00001
%
% FORMAT files_out = spmup_hrf_boost
%        files_out = spmup_hrf_boost(SPM_location,shift,columns_of_interest)
%
% INPUT SPM_location is the full name of the SPM.mat
%       shift is time shift allowed around the hrf peak - default is 2.5
%             if 0 use, no masking is applied (ie we likely fit noise)
%       columns_of_interest is a vector for the columns of the hrf in X that
%       need boosting (allows skipping non interesting betas) - default = []
%
% OUTPUT files_out list files created in /hrf_boost
%        beta_XXXX corresponding to the boosted hrf
%        con_XXXX combined boosted hrf regressors
%
% Cyril Pernet April 2014
% Updated August 2015 to handle .nii + batch input
% Updated October 2015 fixing data handling + make output for batch
% February 2016 added the columns_of_interest and optimized code for speed
% -----------------------------------------------------------------------
% Copyright (C) spmup team 2016

disp('Boosting parameter estimates')
disp('Ref: Pernet CR (2014) Misconceptions in the use of the General Linear')
disp('Model applied to functional MRI: a tutorial for junior neuro-imagers.')
disp('Front. Neurosci. 8:1. doi: 10.3389/fnins.2014.00001')

newbeta  = '';
newcon   = '';
columns_of_interest = [];
defaults = spm_get_defaults;
img_ext  = defaults.images.format;

current = pwd;
shift = 2; % default is 4 to 6 sec
if nargin == 0
    [t,sts] = spm_select(1,'mat','Select 1st level SPM file');
    if sts == 0
        return
    end
elseif nargin == 1
    t = varargin{1};
elseif nargin == 2
    t = varargin{1};
    shift = varargin{2};
elseif nargin == 3
    t = varargin{1};
    shift = varargin{2};
    columns_of_interest = varargin{3};
end

if iscell(t); t = cell2mat(t); end
if iscell(shift); shift = cell2mat(t); end
if isstruct(shift); shift = getfield(shift,cell2mat(fieldnames(shift))); end

% check the design
% ------------------
load(t);

if strcmp(SPM.xBF.name,'hrf (with time and dispersion derivatives)')
    type = 3;
elseif strcmp(SPM.xBF.name,'hrf (with time derivative)')
    type = 2;
else
    error('spm_hrf_boost only works for designs with 1st and/or 2nd derivatives')
end

[path,file,ext] = fileparts(t);
cd(path);
mkdir('hrf_boost')

% using a reference trial, estimate the range of possible delays
% use this to constrains the hrf boost, since we don't want to boost
% the hrf with a time derivative that fits some artefacts -
% take as possible range +/- the shift parameter
% --------------------------------------------------------------
disp('--------------------------------')
disp('  Estimating acceptable delays  ')
disp('--------------------------------')

% parameters of the design
xBF.dt     = SPM.xBF.dt;
xBF.name   = SPM.xBF.name;
xBF.length = SPM.xBF.length;
xBF.order  = SPM.xBF.order;
xBF        = spm_get_bf(xBF);

% orthogonalise derivative(s) and normalize
xBF.bf =  spm_orth(xBF.bf);
SS     = xBF.bf'*xBF.bf;
for f=1:size(xBF.bf,2)
    xBF.bf(:,f) = xBF.bf(:,f)./sqrt(SS(f,f));
end

% the range of time covered by this model is
BFtime     = 0:xBF.dt:(length(xBF.bf)-1)*xBF.dt;
PeakTime   = BFtime(find(xBF.bf(:,1) == max(xBF.bf(:,1))));
PeakDelays = [PeakTime-shift PeakTime+shift];
if length(unique(PeakDelays)) == 1
    fprintf('skipping time to peak estimation - delay = 0, this is risky \n')
else
    fprintf('parameters will be boosted for estimated responses peaking between %g and %g sec \n',PeakDelays(1),PeakDelays(2));
end

% boost parameters
% ------------------

disp('--------------------------------')
disp('  Re-estimating hrf magnitude   ')
disp('--------------------------------')

index = 1; i = 1;
for s=1:size(SPM.Sess,2) % for each session
    for n = 1:length(SPM.Sess(s).U) % for each condition
        for c = 1:length(SPM.Sess(s).U(n).name) % for each regressor in this condition
            if sum(index == columns_of_interest) == 1 || isempty(columns_of_interest)
                
                clear im
                % name
                if index < 9
                    name = sprintf('boost_beta_000%g',index);
                    im(1,:) = [pwd filesep sprintf('beta_000%g.%s',index,img_ext)];
                    im(2,:) = [pwd filesep sprintf('beta_000%g.%s',index+1,img_ext)];
                    if type == 3 && index == 8
                        im(3,:) = [pwd filesep sprintf('beta_00%g.%s',index+2,img_ext)];
                    elseif type == 3 && index < 8
                        im(3,:) = [pwd filesep sprintf('beta_000%g.%s',index+2,img_ext)];
                    end
                elseif index == 9
                    name = sprintf('boost_beta_000%g',index);
                    im(1,:) = [pwd filesep sprintf('beta_000%g.%s',index,img_ext)];
                    im(2,:) = [pwd filesep sprintf('beta_00%g.%s',index+1,img_ext)];
                    if type == 3
                        im(3,:) = [pwd filesep sprintf('beta_00%g.%s',index+2,img_ext)];
                    end
                elseif (index>=10) && (index<99)
                    name = sprintf('boost_beta_00%g',index);
                    im(1,:) = [pwd filesep sprintf('beta_00%g.%s',index,img_ext)];
                    im(2,:) = [pwd filesep sprintf('beta_00%g.%s',index+1,img_ext)];
                    if type == 3 && index == 98
                        im(3,:) = [pwd filesep sprintf('beta_0%g.%s',index+2,img_ext)];
                    elseif type == 3 && index < 98
                        im(3,:) = [pwd filesep sprintf('beta_00%g.%s',index+2,img_ext)];
                    end
                elseif index == 99
                    name = sprintf('boost_beta_00%g',index);
                    im(1,:) = [pwd filesep sprintf('beta_00%g.%s',index,img_ext)];
                    im(2,:) = [pwd filesep sprintf('beta_0%g.%s',index+1,img_ext)];
                    if type == 3
                        im(3,:) = [pwd filesep sprintf('beta_0%g.%s',index+2,img_ext)];
                    end
                elseif (index>=100) && (index<999)
                    name = sprintf('boost_beta_0%g',index);
                    im(1,:) = [pwd filesep sprintf('beta_0%g.%s',index,img_ext)];
                    im(2,:) = [pwd filesep sprintf('beta_0%g.%s',index+1,img_ext)];
                    if type == 3 && index == 998
                        im(3,:) = [pwd filesep sprintf('beta_%g.%s',index+2,img_ext)];
                    elseif type == 3 && index < 998
                        im(3,:) = [pwd filesep sprintf('beta_0%g.%s',index+2,img_ext)];
                    end
                elseif index == 999
                    name = sprintf('boost_beta_0%g',index);
                    im(1,:) = [pwd filesep sprintf('beta_0%g.%s',index,img_ext)];
                    im(2,:) = [pwd filesep sprintf('beta_%g.%s',index+1,img_ext)];
                    if type == 3
                        im(3,:) = [pwd filesep sprintf('beta_%g.%s',index+2,img_ext)];
                    end
                elseif index>=1000
                    name = sprintf('boost_beta_%g',index);
                    im(1,:) = [pwd filesep sprintf('beta_%g.%s',index,img_ext)];
                    im(2,:) = [pwd filesep sprintf('beta_%g.%s',index+1,img_ext)];
                    if type == 3
                        im(3,:) = [pwd filesep sprintf('beta_%g.%s',index+2,img_ext)];
                    end
                end
                
                
                % inform user
                if type == 2
                    fprintf('Combining \n %s \n %s \n',im(1,:),im(2,:));
                elseif type == 3
                    fprintf('Combining \n %s \n %s \n %s \n',im(1,:),im(2,:),im(3,:));
                end
                
                % Create a mask image for voxels outside the PeakDelays
                V = spm_vol(im);
                sts = spm_check_orientations(V);
                images = spm_read_vols(V);
                
                if length(unique(PeakDelays)) == 1
                    Mask = ~isnan(squeeze(images(:,:,:,1)));
                else
                    T2P = NaN(V(1).dim); % time to peak image
                    voxels = find(~isnan(squeeze(images(:,:,:,1))));
                    [x,y,z] = ind2sub(size(squeeze(images(:,:,:,1))),voxels);
                    time = nan(1,length(voxels));
                    BF = SPM.xBF.bf;
                       parfor v=1:length(voxels)
                        if images(x(v),y(v),z(v),1) > 0 % if positive beta take the max
                            if type == 2
                                time(v) = BFtime(find((BF*[images(x(v),y(v),z(v),1) images(x(v),y(v),z(v),2)]')==max((BF*[images(x(v),y(v),z(v),1) images(x(v),y(v),z(v),2)]'))));
                            elseif type == 3
                                time(v) = BFtime(find((BF*[images(x(v),y(v),z(v),1) images(x(v),y(v),z(v),2) images(x(v),y(v),z(v),3)]')==max((BF*[images(x(v),y(v),z(v),1) images(x(v),y(v),z(v),2) images(x(v),y(v),z(v),3)]'))));
                            end
                            
                        else % take the min
                            if type == 2
                                time(v) = BFtime(find((BF*[images(x(v),y(v),z(v),1) images(x(v),y(v),z(v),2)]')==min((BF*[images(x(v),y(v),z(v),1) images(x(v),y(v),z(v),2)]'))));
                            elseif type == 3
                                time(v) = BFtime(find((BF*[images(x(v),y(v),z(v),1) images(x(v),y(v),z(v),2) images(x(v),y(v),z(v),3)]')==min((BF*[images(x(v),y(v),z(v),1) images(x(v),y(v),z(v),2) images(x(v),y(v),z(v),3)]'))));
                            end
                        end
                       end
                    
                    for v=1:length(voxels); T2P(voxels(v)) = time(v); end
                    
                    newV = spm_create_vol(V(1));
                    [path,file,ext] = fileparts(V(1).fname);
                    newV.fname = [path filesep 'hrf_boost' filesep 'T2P' num2str(index) ext];
                    newV.descrip = sprintf('time to peak image of %s',V(1).descrip);
                    spm_write_vol(newV,T2P);
                    Mask = logical((T2P>=PeakDelays(1)).*(T2P<=PeakDelays(2)));
                end
                
                % boost the hrf values
                X = SPM.xX.X(:,[index index+1 index+2]);
                SumOfSquares = diag(X'*X); % these are the weights
                if sum(SumOfSquares - 1) > 0.0001 % also check X is normalized
                    X(:,1) = X(:,1)./norm(X(:,1));
                    X(:,2) = X(:,2)./norm(X(:,2));
                    X(:,3) = X(:,3)./norm(X(:,3));
                    SumOfSquares = diag(X'*X);
                end
                
                % Create and execute the spm image calculation.
                if type == 2
                    magnitude = sqrt(((squeeze(images(:,:,:,1)).^2).*SumOfSquares(1)) + ...
                        ((squeeze(images(:,:,:,2)).^2).*SumOfSquares(2))); % like Steffener
                elseif type == 3
                    magnitude = sqrt(((squeeze(images(:,:,:,1)).^2).*SumOfSquares(1)) + ...
                        ((squeeze(images(:,:,:,2)).^2).*SumOfSquares(2)) + ...
                        ((squeeze(images(:,:,:,3)).^2).*SumOfSquares(3)));
                end
                sign = (squeeze(images(:,:,:,1))) ./ abs((squeeze(images(:,:,:,1)))); % get sign like Calhoun
                boosted_values = magnitude .* sign;
                
                % Only change values within Mask
                boosted_hrf = squeeze(images(:,:,:,1)); % beta hrf
                boosted_hrf(Mask) = boosted_values(Mask); % only update inside the mask
                
                % cleanup
                try
                    Mask = spm_read_vols(spm_vol([SPM.swd filesep 'mask.' img_ext]));
                catch diff_path
                   Mask = spm_read_vols(spm_vol([pwd filesep 'mask.' img_ext]));
                end
                boosted_hrf = boosted_hrf.*Mask;
                
                % save file
                newV = spm_create_vol(V(1));
                [path,file,ext] = fileparts(V(1).fname);
                newV.fname = [path filesep 'hrf_boost' filesep name ext];
                newbeta(i,:) = [path filesep 'hrf_boost' filesep name ext];
                newV.descrip = sprintf('Boosted version of %s',V(1).descrip);
                spm_write_vol(newV,boosted_hrf);
                fprintf('boost_%s done \n',file);
                
                % update index for the next hrf
                hrf_indices(i) = index; i=i+1;
                if type == 2
                    index = index+2;
                elseif type == 3
                    index = index+3;
                end
                clear name im V images X magnitude sign boosted_hrf newV
            end
        end
    end
    % after a given session, we might have to move further if motion etc
    % are added to the design
    index = SPM.Sess(s).col(end)+1;
end

% also update constrasts
% ----------------------

if isfield(SPM,'xCon')
    
    disp(' ')
    disp('--------------------------------')
    disp('    Re-estimating contrasts     ')
    disp('--------------------------------')
    cd hrf_boost/
    
    index = 1;
    for c = 1:length(SPM.xCon) % for each contrast
        columns = find(SPM.xCon(c).c); % check columns
        test = intersect(columns,hrf_indices);
        if length(test) == length(columns) % constrast involving combination of hrf only
            % load boosted parameters
            for i=1:length(columns)
                if columns(i) < 10
                    im(i,:) = [pwd filesep sprintf('boost_beta_000%g.%s',columns(i),img_ext)];
                elseif (columns(i) >= 10) && (columns(i) < 100)
                    im(i,:) = [pwd filesep sprintf('boost_beta_00%g.%s',columns(i),img_ext)];
                elseif (columns(i) >= 100) && (columns(i) < 1000)
                    im(i,:) = [pwd filesep sprintf('boost_beta_0%g.%s',columns(i),img_ext)];
                elseif columns(i) >= 1000
                    im(i,:) = [pwd filesep sprintf('boost_beta_%g.%s',columns(i),img_ext)];
                end
            end
            
            try
                V = spm_vol(im);
                sts = spm_check_orientations(V);
                images = spm_read_vols(V);
                
                % compute the constrast
                fprintf('Contrasting \n')
                boosted_con = zeros(size(images,1),size(images,2),size(images,3));
                for i=1:length(columns)
                    fprintf('%s \n',im(i,:));
                    boosted_con = boosted_con + (SPM.xCon(c).c(columns(i)).*squeeze(images(:,:,:,i)));
                end
                
                % save file
                newV = spm_create_vol(V(1));
                [path,file,ext] = fileparts(SPM.xCon(c).Vcon.fname);
                newV.fname = [pwd filesep 'boost_' file ext];
                newcon(index,:) = [pwd filesep 'boost_' file ext];
                newV.descrip = sprintf('Boosted version of %s',SPM.xCon(c).Vcon.descrip);
                spm_write_vol(newV,boosted_con); index = index + 1;
                fprintf('boost_%s done \n',SPM.xCon(c).Vcon.fname);
            
            catch spm_vol_err
                if sum(strfind(spm_vol_err,'exist')) ~= 0
                    fprintf('%s skipped, betas were not boosted \n',SPM.xCon(c).Vcon.fname);
                end
            end

        end
    end
end

files_out{1} = newbeta;
files_out{2} = newcon;

cd(current)
disp('------------------')
disp('spm_hrf_boost done')
