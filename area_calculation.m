
%Code for cluster area calculation
% clear all;
close all;
 
 
no_of_clusters=size(clusters,1);
for i=1:no_of_clusters
    data=[];
 data(:,:) = clustersXY{i,1};
 
 
 
% Calculate the eigenvectors and eigenvalues
covariance = cov(data);
[eigenvec, eigenval ] = eig(covariance);
 
% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);
 
% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));
 
% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2))
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1))
    smallest_eigenvec = eigenvec(1,:);
end
 
% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));
 
% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end
 
% Get the coordinates of the data mean
avg = mean(data);
 
% Get the 95% - 2.4474   90- 2.14706    80-1.7944    confidence interval error ellipse
chisquare_val = 2.4474;
theta_grid = linspace(0,2*pi);
phi = angle;
X0=avg(1);
Y0=avg(2);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);
 
%area
area(i,1)=phi*a*b*0.06*0.06;%0.06 is the pixel size in micro meters
number_mol(i,1)=size(data,1);
density(i,1)= number_mol(i,1)/area(i,1);
end

figure;  
bar(sort(area));
ax = gca; 
ax.FontSize = 20;
ax.FontWeight = 'bold';
xlabel('Clusters');
ylabel('Area (\mum^{2})');
  
figure;  
bar(sort(number_mol));
ax = gca; 
ax.FontSize =20;
ax.FontWeight = 'bold';
xlabel('Clusters');
ylabel('Molecules');

 
figure;  
bar(sort(density));
ax = gca; 
ax.FontSize = 20;
ax.FontWeight = 'bold';
xlabel('Clusters');
ylabel('Density');
 

 
 
 
 
%%%%%mean
 
mean_area = sum(area)/no_of_clusters;
mean_mol=sum(number_mol)/no_of_clusters;
mean_density=sum(density)/no_of_clusters;
 
mini_mol=min(number_mol);
maxi_mol=max(number_mol);
