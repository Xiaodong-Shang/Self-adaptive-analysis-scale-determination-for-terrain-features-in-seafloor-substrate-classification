
%% The following is Freedman's method, from paper Multi-Scale Measures of Rugosity, Slope and Aspect from Benthic Stereo Image Reconstructions
function [EastD, NorthD, rgsty, slope, aspect, normal, rangez, sdevz, rgstyXY, meanxyz, nnodes, concavity, meandevz, trinormals, triareas, badinds] = trisurfterrainfeats(tri, xyz, window, wincntrs)
    
    % Calculate areas and normals for every triangle in the tri matrix
    [triareas, trinormals] = trinormarea(tri, xyz);
    
%     spmat = sparse(tri(:,[1 1 2 2 3 3]),tri(:,[2 3 1 3 1 2]),1,numel(xyz(:,1)),numel(xyz(:,1))); % sparse matrix for dsearch
    badcount = 0;
    badinds = [];
    nclrchar = 0; % keep track of printed chars (pretty display purposes).
    
    % output features for entire surface
    if nargin == 2
        % position of centres of triangles
        trixyz1 = xyz(tri(:,1),:); 
        trixyz2 = xyz(tri(:,2),:); 
        trixyz3 = xyz(tri(:,3),:);
        trixyz = [mean([trixyz1(:,1) , trixyz2(:,1) , trixyz3(:,1)],2) , ...
            mean([trixyz1(:,2) , trixyz2(:,2) , trixyz3(:,2)],2) , ...
            mean([trixyz1(:,3) , trixyz2(:,3) , trixyz3(:,3)],2)];

        [rgsty, slope, aspect, normal, rangez, sdevz, rgstyXY, meanxyz, nnodes, concavity] = get_features(xyz,trinormals,triareas,trixyz);
        
        if nargout > 10,
            %zind = dsearch(xyz(:,1),xyz(:,2),tri,meanxyz(1),meanxyz(2),spmat);
            [minzdist zind] = min((xyz(:,1)-meanxyz(1)).^2+(xyz(:,2)-meanxyz(2)).^2);
            meandevz = xyz(zind,3) - meanxyz(3);
        end

%         [normal, slope, aspect] = pcaxyz(xyz,k);
%         triProjAreas = abs(triareas.*dot(trinormals, repmat(normal , size(trinormals,1), 1) , 2));
%         rgsty = sum(triareas)/sum(triProjAreas);
%         
%         % additional outputs
%         rangez = range(xyz(:,3)); 
%         sdevz = std(xyz(:,3)); 
%         triProjAreasXY = abs(triareas.*dot(trinormals, repmat([0 0 -1] , size(trinormals,1), 1) , 2));
%         rgstyXY = sum(triareas)/sum(triProjAreasXY);
%         meanxyz = mean(xyz,1); 
%         nnodes = size(xyz,1); 
        
    % iterate over wincntrs calculating features over window
    elseif nargin >= 3                  
        % if wincntrs not specified, compute for every node in the mesh
        if nargin == 3, wincntrs = xyz(:,[1,2]); end
        
        % Some variable used for the loop
        n = size(wincntrs,1);       % number of elements in wincntrs
        wbtic = round(n/100);    % the number of iterations between updating the progress display
        
        % Preallocate some variables for speed
        rgsty = nan(n,1);     % rugosity variable
        slope = nan(n,1);     % slope variable
        aspect = nan(n,1);     % aspect variable
        rangez = nan(n,1);     % range variable
        sdevz = nan(n,1);     % std dev variable
        rgstyXY = nan(n,1);     % xyplane rugosity variable
        normal = nan(n,3);     % normal variable
        meanxyz = nan(n,3);     % coordinates of mean of points within window
        nnodes = nan(n,1);     % number of xyz points included in window
        concavity = nan(n,1);     % number of concavity points included in window
        meandevz = nan(n,1);     % number of relative zdev points included in window
                
        % Iterate through every point in wincntrs
        for i = 1:n
            try
                % Generate mask for points within the specified search window
                [trimsk, xyzmsk] = trisurfcropmask(tri, xyz, window, wincntrs(i,:));   % obtain mask for tri matrix
                %tri_subset = tri(trimsk,:);
                xyz_subset = xyz(xyzmsk, :);

                % position of centres of triangles
                trixyz1 = xyz(tri(trimsk,1),:); trixyz2 = xyz(tri(trimsk,2),:); trixyz3 = xyz(tri(trimsk,3),:);
                trixyz = [mean([trixyz1(:,1) , trixyz2(:,1) , trixyz3(:,1)],2) , ...
                    mean([trixyz1(:,2) , trixyz2(:,2) , trixyz3(:,2)],2) , ...
                    mean([trixyz1(:,3) , trixyz2(:,3) , trixyz3(:,3)],2)];

                [rgsty(i), slope(i), aspect(i), normal(i,:), rangez(i), sdevz(i), rgstyXY(i), meanxyz(i,:), nnodes(i), concavity(i)] = get_features(xyz_subset,trinormals(trimsk,:),triareas(trimsk),trixyz,wincntrs(i,:));

                if nargout > 10,
                    %zind = dsearch(xyz(:,1),xyz(:,2),tri,meanxyz(i,1),meanxyz(i,2),spmat);
                    [minzdist zind] = min((xyz(:,1)-wincntrs(i,1)).^2+(xyz(:,2)-wincntrs(i,2)).^2);
                    meandevz(i) = xyz(zind,3) - meanxyz(i,3);
                end
            catch ME
                badcount = badcount+1;
                badinds(badcount) = i;
            end

            
            if ~mod(i,wbtic)
                [nclrchar] = clrfprintf(nclrchar , '%d/%d (%.0f%%)',i,n,i/n*100);
            end
        end
        clrfprintf(nclrchar , '');
        if (badcount>0)
            warning([num2str(badcount),' iterations were dodgy.'])
        end
    end
    EastD=sin(aspect);
    NorthD=cos(aspect);
end


function [rgsty, slope, aspect, normal, rangez, sdevz, rgstyXY, meanxyz, nnodes, concavity] = get_features(xyz,trinormals,triareas,trixyz,cntr)
    % unit vector pointing in positive down z-direction
    %k = [0 0 -1]; 
    k = [0 0 -1];
    
    % Do PCA to fit plane    
    [normal, slope, aspect] = pcaxyz(xyz,k);
    
    % Simple outputs
    rangez = range(xyz(:,3)); 
    sdevz = std(xyz(:,3)); 
    nnodes = size(xyz,1); 
    meanxyz = mean(xyz,1); 

    % Calculate rugosity
    triProjAreas = abs(triareas.*dot(trinormals, repmat(normal , size(trinormals,1), 1) , 2));
    rgsty = sum(triareas)/sum(triProjAreas);
    
    % Calculate horizontal plane rugosity
    triProjAreasXY = abs(triareas.*dot(trinormals, repmat([0 0 -1] , size(trinormals,1), 1) , 2));
    rgstyXY = sum(triareas)/sum(triProjAreasXY);

    % Calculate concavity
    if nargin < 5
        cntr = meanxyz(1:2);
    else
        cntr = cntr(1:2);
    end
    vect = repmat(cntr , size(trinormals,1),1) - [trixyz(:,1) trixyz(:,2)];
    normvect = vect./repmat((vect(:,1).^2+vect(:,2).^2).^0.5 , 1 ,2);
%     concavity = mean(dot(trinormals(:,[1, 2]),normvect , 2))./size(trinormals,1);
    concavity = mean(dot(trinormals(:,[1, 2]),normvect , 2));
end


% TRISURFCROPMASK Mask of verticies in a triangular surface in x-y search range.
%   MSK = TRISURFCROPMASK(TRI, XYZ, XRANGE, YRANGE) provides a mask, MSK, for 
%   TRI that contains verticies that lie within the search bounds specified
%   by XRANGE and YRANGE. XRANGE and YRANGE are both 2 element arrays
%   containing the min and max values for x and y, respectively. 
%   
%   Example:
%   [x,y,z] = peaks(40); 
%   tri = delaunay(x,y);
%   xyz = [reshape(x, numel(x), 1) , reshape(y, numel(y), 1) , reshape(z, numel(z), 1)];
%   mask = trisurfcropmask(tri, xyz, [-1,2], [-1,2]);
%   figure, trisurf(tri,xyz(:,1),xyz(:,2),xyz(:,3)), hold on
%   trisurf(tri(mask,:),xyz(:,1),xyz(:,2),xyz(:,3)+10)
%   
%   TODO: update documentation for additional output arg, xyzmsk

%   Author:     Ariell Friedman (a.friedman@acfr.usyd.edu.au)
%   Created:    15/09/2009
%   Updated:    -

function [trimsk, xyzmsk] = trisurfcropmask(tri, xyz, window, loc)
    if length(window) == 1
        % Use axis aligned square window (quicker)
        xrange = [loc(1)-window/2 , loc(1)+window/2];
        yrange = [loc(2)-window/2 , loc(2)+window/2];
        xyzmsk = (xyz(:,1)>=xrange(1) & xyz(:,1)<=xrange(end) & ...  % x-mask
            xyz(:,2)>=yrange(1) & xyz(:,2)<=yrange(end));            % y-mask
    elseif length(window) == 2
        % Use axis aligned rectangular window (quicker)
        xrange = [loc(1)-window(1)/2 , loc(1)+window(1)/2];
        yrange = [loc(2)-window(2)/2 , loc(2)+window(2)/2];
        xyzmsk = (xyz(:,1)>=xrange(1) & xyz(:,1)<=xrange(end) & ...  % x-mask
            xyz(:,2)>=yrange(1) & xyz(:,2)<=yrange(end));            % y-mask
    elseif size(window,1) > 3 && size(window,2)==2
        % Use polygon shaped window (slower)
        % TODO: apply rotation on window based on auv heading
        % TODO: apply changing shape of polygon based on position?
        srchwin = [window(:,1)+loc(1) , window(:,2)+loc(2)];
        %xyzmsktransp = inpoly(xyz(:,[1,2])',srchwin');
        %xyzmsk = xyzmsktransp';
        xyzmsk = inpolygon(xyz(:,1),xyz(:,2),srchwin(:,1),srchwin(:,2));
    end
    
    indsall = 1:size(xyzmsk,1);
    inds = indsall(xyzmsk);
    
    trimskxyz = ismember(tri, inds);
    trimsk = trimskxyz(:,1) & trimskxyz(:,2) & trimskxyz(:,3);
end




function [normal, slope, aspect] = pcaxyz(xyz,dirvect)
    coeff = pca(xyz,'economy', false);
    normal = coeff(:,3)';
        
    if dot(normal, dirvect) < 0, normal=-1*normal; end    % enforce that normal has an upward facing component
    
    %slope = acos(dot(normal, dirvect)/norm(normal)*norm(dirvect));
    %slope = abs(atan(coeff(3,1)/(coeff(1,1)^2+coeff(2,1)^2)))
    slope = abs(acos(dot(normal, dirvect)));
    aspect = atan2(normal(2),normal(1));
end


% TRIAREA3 Area of the triangles in a triangular surface.
%   A = TRIAREA3(TRI, XYZ) Calculates the surface area of the triangular
%   surfaces defined by TRI and XYZ.
%   XYZ is an Mx3 matrix that contains x,y and z coordinates to every
%   vertex in the surface. TRI is an Nx3 matrix where each row contains 
%   indecies to the XYZ vertex matrix to define a single triangular face.
%   The output A is an Nx1 containing the area of every triangle in TRI.
%   
%   The areas are calculated to be half the magnitude of the cross product 
%   of the vectors representing two adjacent sides of the triangle. The 
%   basis for the calculation is as follows: Define the triangle by points 
%   PQR. The area of the parallelogram with sides PQ and PR is equal to the
%   magnitude of the cross product of vectors representing two adjacent 
%   sides. The area of the triangle is then half of this.
%   
%   To obtain the area of the full surface, sum the result. For example:
%       SA = sum(triarea3(tri, xyz));
%
%   TODO: Update documentation
   
%   Author:     Ariell Friedman (a.friedman@acfr.usyd.edu.au)
%   Created:    15/09/2009
%   Updated:    -

function [A, N] = trinormarea(tri, xyz)
    P = xyz(tri(:,1),:);
    Q = xyz(tri(:,2),:);
    R = xyz(tri(:,3),:);

    PQ = Q - P;
    PR = R - P;
    
    N = cross(PQ,PR);               % normals of triangles (X-prod)
    magN = sum(N.^2 , 2).^0.5;      % magnitude of X-prod
    
    A = 0.5 .* magN;                % areas of triangles
    N = N./repmat(magN,1,size(N,2));    % unit vector normals
    
%     k = repmat([0 0 1],size(N,1),1); % unit vector pointing in upward z-direction
%     msk = dot(N, k, 2) < 0;
%     msk = msk*-1;
%     msk(msk ==0) = ones(size(msk(msk ==0)));
%     N = repmat(msk, 1, size(N,2)).*N;
end


% Screen print function. Same as fprintf with the additional first input
% parameter to specify how many characters to delete before printing.
% Useful for neat, self deleting displays.
function [nclrchar] = clrfprintf(nclrchar , varargin)
    fprintf(repmat('\b',1,nclrchar));
    nclrchar = fprintf(varargin{:});
end