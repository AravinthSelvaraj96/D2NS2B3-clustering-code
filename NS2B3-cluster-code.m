%Code for cluster analysis. 8march 2023 pagfp ns2b3 
 

 
%XY= importdata('/Users/aravinth/Documents/MATLAB/final_p1/90ms_xy_new.xlsx'); % read X and Y

XY(:,1)=G_arr(:,5);%xy coordinates of the molecules
XY(:,2)=G_arr(:,6);
 
clearvars -except XY G_arr
figure,plot(XY(:,1), XY(:,2), '.k');


Distance = pdist(XY,'euclidean'); % compute pairwise distances between all points
D_transpose = Distance';
kmax = size(XY,1);
 
%% parameters
maxdist = 1; % distance interms of pixels 
maxdist2 = 1; %distance interms of pixels 
minClusterSize =50; 
 
% Clustering
[Ds,IX] = sort(Distance,'ascend'); % sorts the distances in ascending order
D_sorted_transpose = Ds';
IX_transpose = IX';
 
% Preallocating the variable pairs for speed (thanks to CARLOS RANZ)
indexPairs = 0;
index = 1;
 
 
for k = 1:kmax
    indexPairs = indexPairs + kmax-k;
    if (k < kmax)
        index = [index; indexPairs+1];
    end
end
pairs = ones(indexPairs, 2); % column 1 = point1, column2 = point2
for k = 2:length(index )
    pairs(index(k-1):index(k)-1,:) = [ (k:kmax)' , (k-1) * ones(kmax-(k-1),1 ) ];
end
pairs_sorted = pairs(IX,:);
 
clusters = {}; % list of clusters (one cell per cluster)
pointsInClusters = []; % keep track of which points are already in a cluster (row 1 = cluster number, row 2 = points in cluster)
inc = 0;
tic% increment to keep track of the cluster number
for k = 1:numel(pairs(:,1))
    pt1 = pairs_sorted(k,1);
    pt2 = pairs_sorted(k,2);
    d = Ds(k); % distance between pt1 and pt2
   
    if d <= maxdist
        if isempty(pointsInClusters) || ...
                ( ~ismember(pt1,pointsInClusters(2,:)) && ...
                ~ismember(pt2,pointsInClusters(2,:)) )
           
            inc = inc + 1; % increment + 1 cluster
           
            currentCluster = [pt1,pt2];
            clusters = [ clusters ; {currentCluster} ]; % create new cluster.
            pointsInClusters = [ pointsInClusters , ...
                [ inc * ones(1,length(currentCluster)) ; currentCluster ] ];
           
        elseif ~ismember(pt1,pointsInClusters(2,:))
            id = pointsInClusters(1,pointsInClusters(2,:)==pt2); % find the number of the previous cluster where the point 'pt2' already exists.
            cdist = sqrt( sum((mean(XY(clusters{id},:),1) - XY(pt1,:)).^2) ); % distance between 'pt1' and cluster centroid.
           
            if isempty(cdist) || cdist <= maxdist % add 'pt1' to cluster if within 'maxdist' distance of centroid.
                clusters{id} = [ clusters{id} , pt1 ]; % add 'pt1' to the same cluster, in which 'pt2' is.
                pointsInClusters = [ pointsInClusters , [ id ; pt1 ] ]; % add 'pt1' to list of 'points in clusters'.
            end
           
        elseif ~ismember(pt2,pointsInClusters(2,:))
            id = pointsInClusters(1,pointsInClusters(2,:)==pt1); % find the number of the previous cluster where the point 'pt1' already exists.
            
            cdist = sqrt( sum((mean(XY(clusters{id},:),1) - XY(pt2,:)).^2) ); % distance between 'pt2' and cluster centroid.
                
           
            if isempty(cdist) || cdist <= maxdist % add 'pt2' to cluster if within 'maxdist' distance of centroid or median point.
                clusters{id} = [ clusters{id} , pt2 ]; % add 'pt2' to the same cluster, in which 'pt1' is.
                pointsInClusters = [ pointsInClusters , [ id ; pt2 ] ]; % add 'pt2' to list of 'points in clusters'.
            end
           
        else % if both points are already in different clusters
%             if exist('mergeflag','var')==1 && strcmpi(mergeflag,'merge') % merge both clusters
                id1 = pointsInClusters(1,pointsInClusters(2,:)==pt1); % find the number of the previous cluster where the point 'pt1' already exists.
                id2 = pointsInClusters(1,pointsInClusters(2,:)==pt2); % find the number of the previous cluster where the point 'pt2' already exists.
                if id1~=id2 % merge clusters if different
                    
                    cdist = sqrt( sum((mean(XY(clusters{id1},:),1) - mean(XY(clusters{id2},:),1)).^2) ); % distance between cluster centroids.
                        
                   
                    if isempty(cdist) || cdist <= maxdist2 % merge clusters if within 'maxdist' distance.
                        clusters{id1} = [ clusters{id1} , clusters{id2} ]; % merge both clusters.
                        clusters{id2} = {}; % empty the second cluster (but do not delete it or the cell indices will be messed up!).
                        pointsInClusters(1,pointsInClusters(1,:)==id2) = id1; % replace id of cluster 2 with id of cluster cluster 1.
                    end
                end
%             else % do not merge clusters
                % do nothing
            end
   
   
       else % distance is greater than 'maxdist'
        break; % exit the 'for' loop
    end
       
end
 
clearvars Distance Ds D_transpose D_sorted_transpose IX IX_transpose pairs pairs_sorted
 

toc;
clusters = clusters(cellfun(@(clusters) numel(clusters) >= minClusterSize, clusters)); % exclude clusters that are smaller than 'minClusterSize'.
% clearvars pointsInClusters; % variable is not needed (and not up-to-date) past this point
% Make list of points that are not in a cluster
pointsNotInClusters = (1:kmax);
pointsNotInClusters = pointsNotInClusters(~ismember(pointsNotInClusters,cell2mat(clusters'))); % select points that are not already in a cluster
if minClusterSize < 2 % convert single points into 'one-point' clusters
    clusters = [ clusters ; num2cell(pointsNotInClusters)' ]; % create new clusters.
    pointsNotInClusters= []; % all points are not in a cluster
end
 
% Print summary values
clusterCount = numel(clusters); % Total number of clusters
sizeOfClusters = cellfun(@(clusters) numel(clusters), clusters);
clusterMinSize = min(sizeOfClusters);
clusterMaxSize = max(sizeOfClusters);
clusterMeanSize = mean(sizeOfClusters);
clusterMedianSize = median(sizeOfClusters);
singlepointCount = numel(pointsNotInClusters);
fprintf(1,'Number of clusters: %lu\n',clusterCount); % display the number of clusters
fprintf(1,'Size of smallest cluster: %lu\n',clusterMinSize); % display the number of clusters
fprintf(1,'Size of largest cluster: %lu\n',clusterMaxSize); % display the number of clusters
fprintf(1,'Mean cluster size: %f\n',clusterMeanSize); % display the number of clusters
fprintf(1,'Median cluster size: %lu\n',clusterMedianSize); % display the number of clusters
fprintf(1,'Number of points that are not part of any cluster: %lu\n',singlepointCount); % display the number of single points
% Write the XY coordinates of every points of every cluster
clustersXY = cell(clusterCount,1);
for k=1:clusterCount
    clustersXY{k,1} = XY(clusters{k,1},:);
end
% Compute the centroid (geometrical mean) of every cluster
clustersCentroids = NaN(clusterCount,2);
for k=1:clusterCount
    clustersCentroids(k,:) = mean(clustersXY{k,1},1);
end
% Compute the geometric median of every cluster
clustersGeoMedians = NaN(clusterCount,2);
for k=1:clusterCount
    clustersGeoMedians(k,:) = median(clustersXY{k,1},1);
end
 
 
% Plot the clusters
cc=hsv(clusterCount); % create different colour codes for every cluster
cc = cc(randperm(clusterCount),:); % randomise the colour codes so that neighbouring clusters don't have too similarly looking colours
%h1 = figure('Name','Clusters');
%figure,plot(XY(:,1), XY(:,2), '.k','MarkerSize',10);
figure,scatter(XY(:,1),XY(:,2),5,'filled','o','CData',[0.8,0.8,0.8]); % plot the original points in light grey
hold on;
 
for k=1:clusterCount
    plot(clustersXY{k,1}(:,1),clustersXY{k,1}(:,2),'.','Color',cc(k,:),'MarkerFaceColor',cc(k,:),'MarkerSize',10);
end
%axis equal;
% Plot the centroid of each cluster
%h2 = figure('Name','Centroids of clusters');
hold on;
%scatter(XY(:,1),XY(:,2),20,'filled','o','CData',[.8,.8,.8]); % plot the original points in light grey
%scatter(clustersCentroids(:,1),clustersCentroids(:,2),40,'filled','o','CData',cc); % plot the centroids
plot(clustersCentroids(:,1),clustersCentroids(:,2),'k+',...
     'MarkerSize',7,'LineWidth',3)
ax=gca;
ax.FontSize =20;
ax.FontWeight='bold';
box on;
