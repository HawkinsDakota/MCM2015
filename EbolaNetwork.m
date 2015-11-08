DT = 1;
DURATION = 52;
R0 = 2/3;
QUARANTINE = [0.8, 1, 1, 1, 1, 1, 1, 1, 1, 1];
OMEGA = 0.6/3;
RHO = 0.4/3;
BETA = [0, 0.144, 0.144, 0.144, 0.144, 0.144, 0.144, 0.144, 0.144, 0.144];
START = 1;
PREVENTATIVE_VAC = 0; % 0 if no preventative vaccinations
time = 1;
N = 10;
adjacency = zeros(N,N);
infectArray = zeros(N,1);
infectDay = [1,0,0,0,0,0,0,0,0,0];
currentCondition = zeros(N,3);
overall = zeros(DURATION + 1,N*3);
continuousData = cell(1,N*2);

%Initialize adjacency matrix
adjacency(1,2) = 1;
adjacency(1,3)=1;
adjacency(1,4) = 1;
adjacency(1,5) = 1;
adjacency(1,6) = 1;

adjacency(2,1) = 1;
adjacency(2,6) = 1;
adjacency(3,1) = 1;
adjacency(3,4) = 1;

adjacency(4,1) = 1;
adjacency(4,3) = 1;
adjacency(4,5) = 1;
adjacency(4,8) = 1;
adjacency(4,9) = 1;
adjacency(4,10) = 1;

adjacency(5,1) = 1;
adjacency(5,4) = 1;

adjacency(6,1) = 1;
adjacency(6,2) = 1;
adjacency(6,7) = 1;

adjacency(7,6) = 1;

adjacency(8,4) = 1;
adjacency(8,9) = 1;

adjacency(9,4) = 1;
adjacency(9,8) = 1;

adjacency(10,4) = 1;
%%%%%%%%%%%%%%%%%%%%%%
    
%Initializing first infected city, arbitrarily n1 set to default
infectArray(START,1) = 1;
currentCondition(START,1) = 4750;
currentCondition(START,2) = 250;

%Initialize non-infected cities
currentCondition(2,1) = 3500;
currentCondition(3,1) = 3500;
currentCondition(4,1) = 5000;
currentCondition(5,1) = 3500;
currentCondition(6,1) = 3500;
currentCondition(7,1) = 1000;
currentCondition(8,1) = 5000;
currentCondition(9,1) = 1000;
currentCondition(10,1) = 3500;

count = 0;
for i = 1:10
    overall(time,1 + count) = currentCondition(i,1);
    overall(time,2 + count) = currentCondition(i,2);
    overall(time,3 + count) = currentCondition(i,3);
    count = count + 3;
end
[t,y] = ode45(@Ebola2,[time time + DT],[currentCondition(START,1); currentCondition(START,2); currentCondition(START,3)],[],R0, QUARANTINE(1),0, OMEGA, RHO);
rows = size(t,1);

currentCondition(START,1) = y(rows,1);
currentCondition(START,2) = y(rows,2);
currentCondition(START,3) = y(rows,3);

for i = 1:2*N
    continousData{1,i} = zeros(1,3);
end
continuousData{1,1} = [t; continuousData{1,1}];
continuousData{1,2} = [y; continuousData{1,2}];

time = time + DT;

while time < DURATION
    for i = 1:N
        if infectArray(i,1) == 1
            randomProb = rand(N,1);
            for j = 1:N
                y_m = currentCondition(i,:);
                y_n = currentCondition(j,:);
                prob = adjacency(i,j)*SpreadProp(y_m,y_n,QUARANTINE(i),QUARANTINE(j));
                boolTest = and(randomProb(j) < prob, prob > 0);
                if and(boolTest, infectArray(j) ~= 1)
                    infectArray(j,1) = 1;
                    infectDay(1,j) = time;
                    percInfected = rand*0.05;
                    currentCondition(j,2) = percInfected*currentCondition(j,1);
                    currentCondition(j,1) = currentCondition(j,1) - currentCondition(j,2);
                end
            end
        end
    end
    for k = 1:N
        if infectArray(k,1) == 1
            neighbors = find(adjacency(k,:));
            [t,y] = ode45(@Ebola2,[time time + DT],[currentCondition(k,1); currentCondition(k,2); currentCondition(k,3)],[],R0, QUARANTINE(k),BETA(k), OMEGA, RHO);
            rows = size(t,1);
            currentCondition(k,1) = y(rows,1);
            currentCondition(k,2) = y(rows,2);
            currentCondition(k,3) = y(rows,3);
            continuousData{1,2*k -1} = [t; continuousData{1,2*k-1}];
            continuousData{1,2*k} = [y; continuousData{1,2*k}];
            if PREVENTATIVE_VAC == 1;
                total = size(neighbors,1);
                for l = neighbors
                    if infectArray(neighbors) ~= 1
                         [t,y] = ode45(@Ebola2,[time time + DT],[currentCondition(l,1); currentCondition(l,2); currentCondition(l,3)],[],R0, QUARANTINE(l),BETA(l), OMEGA, RHO);
                         rows = size(t,1);
                         currentCondition(k,1) = y(rows,1);
                         currentCondition(k,2) = y(rows,2);
                         currentCondition(k,3) = y(rows,3);
                         continuousData{1,2*l -1} = [t; continuousData{1,2*l-1}];
                         continuousData{1,2*l} = [y; continuousData{1,2*l}];
                    end
                end
            end
                    
        end
    end
    time = time + DT;
    count = 0;
    for i = 1:N
        overall(time,1 + count) = currentCondition(i,1);
        overall(time,2 + count) = currentCondition(i,2);
        overall(time,3 + count) = currentCondition(i,3);
        count = count + 3;
    end
end