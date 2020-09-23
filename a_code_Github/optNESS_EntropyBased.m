function [DepFea opt_nn_size] = optNESS_EntropyBased(XYZ,k_min,k_max,delta_k,n1,n2)
% This is our tool optNESS for deriving optimal 3D neighborhoods via
% eigenentropy-based scale selection.
%   Firstly, we used the  eigenentropy-based scale selection method according to Weinmann's paper.
%   Secondly, we calculated the bathymetric features according to Evan's, Wilson's and Friedman's method by using the selected points to build a surface.
%   Please cite Weinmann's paper, Evan's paper, Wilson's paper,Friedman's paper and my paper if you used this function.
%

% get point IDs
point_ID_max = n2-n1+1;
k = k_min:delta_k:k_max;
k=k.*k;


% get local neighborhoods consisting of k neighbors
data_pts = XYZ(n1:n2,1:3);
k_plus_1 = max(k);
num_k = length(k);

[idx,~] = knnsearch(XYZ(:,1:2),data_pts(:,1:2),'Distance','euclidean','NSMethod','kdtree','K',k_plus_1);

nNum=length(data_pts);

% do some initialization stuff for incredible speed improvement
% Shannon_entropy = zeros(point_ID_max,num_k);
opt_nn_size = zeros(point_ID_max,1);
DepFea= zeros(point_ID_max,18);
% calculate Shannon entropy
for j=1:nNum
    Shannon_entropy_real = zeros(1,num_k);  
    for j2=1:num_k

        % select neighboring points
        P = XYZ(idx(j,1:k(j2)),:);          % the point and its k neighbors ...
%         P = data_pts(idx0(1:k(j2)),:);  
        [m,~] = size(P);

        % calculate covariance matrix C
        % 1.) Standard Matlab code (quite slow for small matrices):
        %        C = cov(P);
        % 2.) Fast Matlab code:
        P = P-ones(m,1)*(sum(P,1)/m);
        C = P.'*P./(m-1);

        % get the eigenvalues of C (sorting is already done by Matlab routine eig)
        % ... and remove negative eigenvalues (NOTE: THESE APPEAR ONLY BECAUSE OF NUMERICAL REASONS AND ARE VERY VERY CLOSE TO 0!)
        % ... and later avoid NaNs resulting for eigenentropy if one EV is 0
        [~, D] = eig(C);

        epsilon_to_add = 1e-8;
        EVs = [D(3,3) D(2,2) D(1,1)];
        if EVs(3) <= 0; EVs(3) = epsilon_to_add;
            if EVs(2) <= 0; EVs(2) = epsilon_to_add;
                if EVs(1) <= 0; EVs(1) = epsilon_to_add; end;
            end;
        end;

        % normalize EVs
        EVs = EVs./sum(EVs(:));

        % derive Shannon entropy based on eigenentropy
        Shannon_entropy_cal = -( EVs(1)*log(EVs(1)) + EVs(2)*log(EVs(2)) + EVs(3)*log(EVs(3)) );
        Shannon_entropy_real(j2) = real(Shannon_entropy_cal);

    end  % j2       

%     Shannon_entropy(j,:) = Shannon_entropy_real;

    % select k with minimal Shannon entropy
    [~,min_entry_of_Shannon_entropy] = min(Shannon_entropy_real(:));
    opt_nn_size(j,1) = k(min_entry_of_Shannon_entropy);
    
    % grid surface
    clear P
    P = XYZ(idx(j,1:opt_nn_size(j)),:);
    % select features
    East=P(:,1);
    North=P(:,2);
    Depth=P(:,3);
    DepthCenter=data_pts(j,:);
    [BathyFeatures2]=DepthAnalyze(East,North,Depth,DepthCenter);
    tri=delaunay(P(:,1),P(:,2));
    [EastD, NorthD, r, s, a, ~, ~, ~, rXY] = trisurfterrainfeats(tri, P);
    BathyFeatures1=[EastD, NorthD, r, s, a,rXY];
    
    DepFea(j,:)=[BathyFeatures1 BathyFeatures2];

end  % j


end  % function

