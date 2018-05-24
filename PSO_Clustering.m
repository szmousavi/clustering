function [ J, intra_dist, inter_dist ] = PSO_Clustering( data, clusters_n, max_t )

data_n = 2400;
particles_n = 20;
d_n = size(data,2);

c0 = .7;
c1 = 1.5;
c2 = 1.5;


d = zeros(particles_n,data_n,clusters_n);
cluster = zeros(particles_n,data_n);
particles = zeros(particles_n, clusters_n, d_n);
v = zeros(particles_n, clusters_n, d_n);
local_best = ones(particles_n, clusters_n, d_n);
local_best_err = 1000*ones(particles_n, 1);

% initialization
for i=1:particles_n
    particles(i,:,:) = data(randi(data_n, 1, clusters_n),:);
end

for t=1:max_t
    err = zeros(1, particles_n);
    for i=1:particles_n
        for j=1:data_n
           for k=1:clusters_n
              d(i,j,k) = sqrt((data(j,1)-particles(i,k,1))^2+(data(j,2)-particles(i,k,2))^2);
           end
           [~, idx] = min(d(i,j,:));
           cluster(i,j) = idx(1);
        end
        
        for j=1:clusters_n
            err(i) = err(i) + sum(d(i,cluster(i,:)==j,j))/size(cluster(i,:)==j,2);
        end
        err(i) = err(i)/clusters_n;
    end
    
    for i=1:particles_n    
        if err(i) < local_best_err(i)
            local_best(i,:,:) = particles(i,:,:); 
            local_best_err(i) = err(i);
        end
    end
    
    [~, idx] = min(local_best_err);
    
    
    for i=1:particles_n
        r1 = rand(1,clusters_n, d_n);
        r2 = rand(1,clusters_n, d_n);
        v(i,:,:) = c0*v(i,:,:) + c1*r1.*(local_best(i,:,:)-particles(i,:,:)) + c2*r2.*(local_best(idx(1),:,:)-particles(i,:,:));
        particles(i,:,:) = particles(i,:,:) + v(i,:,:);
    end
end

global_best = local_best(idx(1),:,:);
J = local_best_err(idx(1))

d = zeros(data_n,clusters_n);
cluster = zeros(1,clusters_n);
intra_dist = 0;
for j=1:data_n
   for k=1:clusters_n
      d(j,k) = sqrt((data(j,1)-global_best(1,k,1))^2+(data(j,2)-global_best(1,k,2))^2);%sum((data(j,:)-particles(i,k,:)).^2)); 
   end
   [dist, idx] = min(d(j,:));
   intra_dist = intra_dist+dist;
   cluster(j) = idx(1);
end
intra_dist = intra_dist/data_n

inter_dist = 0;
for i=1:clusters_n
    for j=1:clusters_n
        inter_dist = inter_dist+sqrt((global_best(1,i,1)-global_best(1,j,1))^2+(global_best(1,i,2)-global_best(1,j,2))^2);
    end
end
inter_dist = inter_dist/clusters_n^2


% figure;
% plot(data(cluster==1,1),data(cluster==1,2),'r.','MarkerSize',12)
% hold on
% plot(data(cluster==2,1),data(cluster==2,2),'b.','MarkerSize',12)
% plot(data(cluster==3,1),data(cluster==3,2),'g.','MarkerSize',12)
% plot(data(cluster==4,1),data(cluster==4,2),'m.','MarkerSize',12)
% plot(local_best(idx(1),:,1),local_best(idx(1),:,2),'kx',...
%      'MarkerSize',15,'LineWidth',3)
% legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Centroids',...
%        'Location','NW')
% title 'Cluster Assignments and Centroids'
% hold off

end