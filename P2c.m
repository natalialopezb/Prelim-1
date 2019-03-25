%% Prelim 1.2c

%% Case 1 - Steady-state

[m,n] = size(s1);
temp_s1 = zeros(m,1);
for i=1:m
    for j=1:n
        temp_s1(i,1) = temp_s1(i,1)+s1(i,j);
    end
end

BigS1 = zeros(6,m/6);
aux = 0; %For me
for j=1:m/6
    for i=1:6
        BigS1(i,j)=temp_s1(i+aux,1);
    end
    aux = aux+6;
end

[U1,S1,V1] = svd(BigS1);
U1 = abs(U1(:,1));

%% Case 2 - Early Inducer

[m,n] = size(s2);
temp_s2 = zeros(m,1);
for i=1:m
    for j=1:n
        temp_s2(i,1) = temp_s2(i,1)+s2(i,j);
    end
end

BigS2 = zeros(6,m/6);
aux = 0; %For me
for j=1:m/6
    for i=1:6
        BigS2(i,j)=temp_s2(i+aux,1);
    end
    aux = aux+6;
end

[U2,S2,V2] = svd(BigS2);
U2 = abs(U2(:,1));

%% Case 3 - Late Inducer

[m,n] = size(s3);
temp_s3 = zeros(m,1);
for i=1:m
    for j=1:n
        temp_s3(i,1) = temp_s3(i,1)+s3(i,j);
    end
end

BigS3 = zeros(6,m/6);
aux = 0; %For me
for j=1:m/6
    for i=1:6
        BigS3(i,j)=temp_s3(i+aux,1);
    end
    aux = aux+6;
end

[U3,S3,V3] = svd(BigS3);
U3 = abs(U3(:,1));