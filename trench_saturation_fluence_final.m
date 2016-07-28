clear
clc

iter = 1; % Number of simulation iterations for averaging 
          %(SET TO AMOUNT OF WORKERS IN PARALLEL POOL)

%------------------ARRAYS FOR DATA STORAGE - DON'T CHANGE------------------
% Arrays for storing intermediate + end depostition results
% wallDepositions4 = zeros(100,iter);
% wallDepositions5 = zeros(100,iter);
% wallDepositions55 = zeros(100,iter);
% wallDepositions6 = zeros(100,iter);
% wallDepositions56 = zeros(100,iter);
% bottomSaturations4 = zeros(iter,1);
% bottomSaturations5 = zeros(iter,1);
% bottomSaturations55 = zeros(iter,1);
% bottomSaturations6 = zeros(iter,1);
% bottomSaturations56 = zeros(iter,1);

%--------------------------GENERAL PARAMETERS------------------------------
writeToTxt = 1;

%----------------REFLECTION-------------
reflectionType = 1; % 0 for specular, 1 for cosinal

%-----------------PROCESS---------------
s = 10^-1; % sticking probability
physTrenchWidth = 50; % trench width (nm) (MUST BE DIVISIBLE BY 100)
AR = 20; % trench aspect ratio (MUST BE INTEGER)
sitesPerNm = 3; % reactive sites per nm (sqrt 8 ~ 3) (MUST BE INTEGER)

numOfParticles = 10^99; % max number of particles to inject

%--------------------------------------------------------------------------

% convert dimensions to number of elements [DON'T CHANGE]
trenchWidth = physTrenchWidth * sitesPerNm;
trenchHeight = AR * physTrenchWidth * sitesPerNm;

% data storage
wallDepositionsEnd = zeros(trenchHeight/100,iter);
bottomSaturationsEnd = zeros(iter,1);
endNumbers = zeros(100,iter);
        
for j = 1:iter
    
    % Array for intermediate particle number data dumps
    endNumbersIter = zeros(100,1);    
      
    % matrices for continuous storing of deposition
    trenchLWall = zeros(trenchHeight,1);
    trenchRWall = zeros(trenchHeight,1);
    trenchBottom = zeros(1,trenchWidth);
    
    % random pools for particle starting positions & angles
    rAnglePool = 2*rand(10^6+1,1)-1;
    rPositionPool = randi([0,trenchWidth],10^6+1,1);    
    
    %--------------------------PROCESS------------------------------
    
    tic
    for i = 1:numOfParticles

        xpos = rPositionPool(mod(i,10^6)+1);
        ypos = trenchHeight;
        
        angle = (pi/2)*rAnglePool(mod(i,10^6)+1);
        vy = -1; 
        vx = tan(angle);
        
        exit = 0;
        reacted = 0;
        collCount = 0;
        
        while exit~=1 && reacted~=1            
            
            if vx < 0;
                % calc left wall intersect point
                leftWallIsct = ypos-(vy/vx)*xpos;
                
                if leftWallIsct > trenchHeight
                    break % particle leaves trench
                end
                
                % reflection off left wall
                if leftWallIsct <= trenchHeight && leftWallIsct > 0
                    % get new position
                    xpos = 0;
                    ypos = leftWallIsct;
                    
                    % get new velocity
                    if reflectionType == 0
                        vx = -vx;
                    elseif reflectionType == 1
                        angle = asin(2*rand-1);
                        vx = 1;
                        vy = tan(angle);
                    end
                    
                    % check reaction
                    if rand < s && trenchLWall(ceil(ypos)) == 0
                        trenchLWall(ceil(ypos))=1;
                        reacted = 1;                        
                    end
                    continue
                end
                
                % reflection off bottom
                if leftWallIsct < 0
                    % get new position
                    xpos = -leftWallIsct/(vy/vx);
                    ypos = 0;
                    
                    % get new velocity
                    if reflectionType == 0
                        vy = -vy;
                    elseif reflectionType == 1
                        angle = asin(2*rand-1);
                        vx = tan(angle);
                        vy = 1;
                    end
                    
                    % check reaction
                    if rand < s && trenchBottom(ceil(xpos)) == 0
                        trenchBottom(ceil(xpos))=1;
                        reacted = 1;
                    end                    
                    continue
                end
                
            else % vx > 0
                
                % calculate right wall intersect
                rightWallIsct = ypos-(vy/vx)*(xpos-trenchWidth);
                
                if rightWallIsct > trenchHeight
                    break % particle leaves trench
                end
                
                % reflection off right wall
                if rightWallIsct <= trenchHeight && rightWallIsct > 0
                    
                    xpos = trenchWidth;
                    ypos = rightWallIsct;
                    
                    if reflectionType == 0
                        vx = -vx;
                    elseif reflectionType == 1
                        angle = asin(2*rand-1);
                        vx = -1;
                        vy = tan(angle);
                    end
                    
                    if rand < s && trenchRWall(ceil(ypos)) == 0
                        trenchRWall(ceil(ypos))=1;
                        reacted = 1;
                    end
                    continue
                end
                
                if rightWallIsct < 0
                    xpos = trenchWidth-rightWallIsct/(vy/vx);
                    ypos = 0;
                    
                    angle = asin(2*rand-1);
                    vx = tan(angle);
                    vy = 1;
                    
                    if rand < s && trenchBottom(ceil(xpos)) == 0
                        trenchBottom(ceil(xpos)) = 1;
                        reacted = 1;
                    end
                    continue
                end
            end
        end
        
        %{
        ---------------------INTERMEDIATE DATA DUMPS-----------------------
        Dump data either after certain number of particles injected, or 
        after certain coverage fraction of the trench surface.
        %}
        
        % PARTICLE MODE
%         if log10(i)==4
%             wallDepositions4(:,j) = plotProfile(trenchHeight,trenchLWall,trenchRWall);
%             bottomSaturations4(j) = nnz(trenchBottom)/1000;
%         elseif log10(i)==5
%             wallDepositions5(:,j) = plotProfile(trenchHeight,trenchLWall,trenchRWall);
%             bottomSaturations5(j) = nnz(trenchBottom)/1000;
%         elseif log10(i)==6
%             wallDepositions6(:,j) = plotProfile(trenchHeight,trenchLWall,trenchRWall);
%             bottomSaturations6(j) = nnz(trenchBottom)/1000;
%         elseif i == 5*10^6
%             wallDepositions56(:,j) = plotProfile(trenchHeight,trenchLWall,trenchRWall);
%             bottomSaturations56(j) = nnz(trenchBottom)/1000;
%         elseif i == 5*10^5
%             wallDepositions55(:,j) = plotProfile(trenchHeight,trenchLWall,trenchRWall);
%             bottomSaturations55(j) = nnz(trenchBottom)/1000;
%         end

        % SATURATION MODE
        saturation = nnz(trenchLWall)/trenchHeight;
        if saturation > 0.003 && mod(saturation,0.01) < 0.003 && ...
                endNumbersIter(round(saturation/0.01)) == 0
            
            endNumbersIter(round(saturation/0.01)) = i;
            disp([num2str(round(saturation/0.01)),'%'])
        end
        
        if endNumbersIter(97)~=0 && ...
                sum(trenchLWall(1:end/4))/sum(trenchLWall(3*(end/4):end)) > 0.95
            
            endNumbersIter(100) = i;
            break
        end
    end
    
    toc % display simulation time
    
    %wallDepositionsEnd(:,j) = plotProfile(trenchHeight,trenchLWall,trenchRWall);
    %bottomSaturationsEnd(j) = nnz(trenchBottom)/trenchWidth;
    endNumbers(:,j) = endNumbersIter;
    
end

%----------------------------AVERAGING RESULTS-----------------------------

% compute average coverages walls
% avgWall4 = mean(wallDepositions4,2);
% avgWall5 = mean(wallDepositions5,2);
% avgWall55 = mean(wallDepositions55,2);
% avgWall6 = mean(wallDepositions6,2);
% avgWall56 = mean(wallDepositions56,2);
avgWallEnd = mean(wallDepositionsEnd,2);
% compute average coverages bottom
% avgBottom4 = mean(bottomSaturations4);
% avgBottom5 = mean(bottomSaturations5);
% avgBottom55 = mean(bottomSaturations55);
% avgBottom6 = mean(bottomSaturations6);
% avgBottom56 = mean(bottomSaturations56);
avgBottomEnd = mean(bottomSaturationsEnd);
% compute average end numbers
avgEndNumbers = mean(endNumbers,2);

%--------------------WRITING END NUMBERS TO TEXT FILE----------------------
if writeToTxt == 1
    % variables for file header
    titleVars = [AR physTrenchWidth log10(s)];
    % generate filename with current parameters
    filename = sprintf('AR %d width %d s 10^%d.txt',titleVars);
    % open file with writing permission
    fileID = fopen(filename,'w');
    % print header to file
    fprintf(fileID,sprintf('AR %d width %d s 10^%d\r\n',titleVars));
    % print data to file
    fprintf(fileID,'\r\n%1d',round(avgEndNumbers).');
    % close file
    fclose(fileID);
end