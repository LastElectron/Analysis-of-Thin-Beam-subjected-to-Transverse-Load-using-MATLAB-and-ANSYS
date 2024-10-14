%%%%%%% Created by Gaurav and Piyush (2023) %%%%%%%
clc;
clear all;
format shortG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Input Values %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx = 96; %Total Length (mm)
Ly = 48; %Total Height (mm)
P = -500; %Point Load (N)
E1 = 25e3; E2 = 50e3; E3 = 100e3; %Young's Modulus (MPa)
E = [E1 E2 E3];
mu = 0.30; %Poisson's Ratio

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Mesh Generation %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Elements
numberElementsX = 48;
numberElementsY = 24;
numberElements = numberElementsX*numberElementsY;

% Element Size
deltaX = Lx/numberElementsX;
deltaY = Ly/numberElementsY;

% Total Nodes
nodesX = numberElementsX+1;
nodesY = numberElementsY+1;
numberNodes = nodesX*nodesY; 

% GDof: global number of degrees of freedom
GDof = 2*numberNodes; 

% nodal coordinates
node = [];
for j = 1:nodesY
    for i = 1:nodesX
        x = (i-1)*deltaX; y = (j-1)*deltaY;
        node = [node; x y];
    end
end
        
% connectivity
element = [];
for j = 1:numberElementsY
    for i = 1:numberElementsX
        i1 = i+(j-1)*nodesX;
        i2 = i1+1;
        i3 = i2+nodesX;
        i4 = i1+nodesX;
        element = [element; i1 i2 i3 i4];
    end
end

nodeCoordinates = node;
elementNodes = element;
xx = nodeCoordinates(:,1);
yy = nodeCoordinates(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Global Stiffness Matrix %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stiffness = zeros(GDof);
thickness =1;

% Gauss point locations & Weights
gaussLocations = [ -0.577350269189626 -0.577350269189626;
                       0.577350269189626 -0.577350269189626;
                       0.577350269189626  0.577350269189626;
                      -0.577350269189626  0.577350269189626];
gaussWeights = [1;1;1;1];

for e = 1:numberElements
    indice = elementNodes(e,:);
    elementDof = [indice indice+numberNodes];
    ndof = length(indice);
    
        % cycle for Gauss point
    for q = 1:size(gaussWeights,1)
        GaussPoint = gaussLocations(q,:);
        xi = GaussPoint(1);
        eta = GaussPoint(2);
        
        % shape functions and derivatives
        shapeFunction = 1/4*[(1-xi)*(1-eta); (1+xi)*(1-eta);
            (1+xi)*(1+eta); (1-xi)*(1+eta)];
        
        naturalDerivatives = 1/4*[
            -(1-eta), -(1-xi); 1-eta,   -(1+xi);
            1+eta ,    1+xi; -(1+eta),  1-xi];
        
        % Jacobian Matrix
        Jacob = nodeCoordinates(indice,:)'*naturalDerivatives;
        invJacobian = inv(Jacob);
        XYderivatives = naturalDerivatives*invJacobian';

        
        % Strain Displacement matrix (B Matrix)
        B = zeros(3,2*ndof);
        B(1,1:ndof)         = XYderivatives(:,1)';
        B(2,ndof+1:2*ndof)  = XYderivatives(:,2)';
        B(3,1:ndof)         = XYderivatives(:,2)';
        B(3,ndof+1:2*ndof)  = XYderivatives(:,1)';
        
        % Material Matrix (D Matrix)
        if e <= numberElements/3
        D = E(1)/(1-mu^2)*[1 mu 0;mu 1 0;0 0 (1-mu)/2];
        elseif (numberElements/3 < e)&&(e <= 2*numberElements/3)
        D = E(2)/(1-mu^2)*[1 mu 0;mu 1 0;0 0 (1-mu)/2];
        elseif (2*numberElements/3 < e)&&(e <= numberElements)
        D = E(3)/(1-mu^2)*[1 mu 0;mu 1 0;0 0 (1-mu)/2];
        end
        
        % Element stiffness matrix at Particular Gauss Point (K Matrix)
        stiffness(elementDof,elementDof) = ...
            stiffness(elementDof,elementDof) + ...
            B'*D*B*thickness*gaussWeights(q)*det(Jacob);
    end
    
 end
    
% Boundary Conditions 
fixedNodeX = find(nodeCoordinates(:,1)==0);  % fixed in XX
fixedNodeY = find(nodeCoordinates(:,1)==0);  % fixed in YY
prescribedDof = [fixedNodeX; fixedNodeY+numberNodes];

% Force vector (Point load applied at Beam End)
force = zeros(GDof,1);
force(2*numberNodes) = P;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Displacement Vector %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

activeDof = setdiff((1:GDof)', prescribedDof);
stiffnessactive = stiffness(activeDof,activeDof);
U = stiffness(activeDof,activeDof)\force(activeDof);
displacements = zeros(GDof,1);
displacements(activeDof) = U;

% Displacement Print
disp('Displacements')
jj = 1:GDof; format("default")
f = [jj; displacements'];
fprintf('node U\n')
fprintf('%3d %12.8f\n',f)
UX = displacements(1:numberNodes);
UY = displacements(numberNodes+1:GDof);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Stress Matrix at Nodes %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodeLocations = [-1 -1;1 -1;1 1;-1 1];
stressExtr = zeros(numberElements,size(nodeLocations,1),3);
strainExtraaa = zeros(numberElements,size(nodeLocations,1),3);

for e = 1:numberElements
    indice3 = elementNodes(e,:);
    elementDof2 = [ indice3 indice3+numberNodes ];
    nn = length(indice3);

    for q = 1:size(gaussWeights,1)
        pt = nodeLocations(q,:);
        xi = pt(1);
        eta = pt(2);

        % shape functions and derivatives
        shapeFunction = 1/4*[(1-xi)*(1-eta); (1+xi)*(1-eta);
            (1+xi)*(1+eta); (1-xi)*(1+eta)];

        naturalDerivatives = 1/4*[
            -(1-eta), -(1-xi); 1-eta,   -(1+xi);
            1+eta ,    1+xi; -(1+eta),  1-xi];

        % Jacobian Matrix
        Jacob = nodeCoordinates(indice3,:)'*naturalDerivatives;
        invJacobian = inv(Jacob);
        XYderivatives = naturalDerivatives*invJacobian';

        % Material Matrix (D Matrix)         
        if e <= numberElements/3
        D = E(1)/(1-mu^2)*[1 mu 0;mu 1 0;0 0 (1-mu)/2];
        elseif (numberElements/3 < e)&&(e <= 2*numberElements/3)
        D = E(2)/(1-mu^2)*[1 mu 0;mu 1 0;0 0 (1-mu)/2];
        elseif (2*numberElements/3 < e)&&(e <= numberElements)
        D = E(3)/(1-mu^2)*[1 mu 0;mu 1 0;0 0 (1-mu)/2];
        end

        %  B matrix
        B = zeros(3,2*nn);
        B(1,1:nn)       = XYderivatives(:,1)';
        B(2,nn+1:2*nn)  = XYderivatives(:,2)';
        B(3,1:nn)       = XYderivatives(:,2)';
        B(3,nn+1:2*nn)  = XYderivatives(:,1)';

        % Strain Vector (at Nodes)
        strainExtr = B*displacements(elementDof2);
        strainExtraaa(e,q,:) = strainExtr;
        % Stress Vector (at Nodes)
        stressExtr(e,q,:) = D*strainExtr;
    end
end

% stress averaging at nodes
stressAvg = zeros(numberNodes,3);
for i = 1:3
    currentStress = stressExtr(:,:,i);
    for n = 1:numberNodes
        idx = find(n==elementNodes);
            stressAvg(n,i) = sum(currentStress(idx))/length(currentStress(idx));
    end
end

% strain averaging at nodes
strainAvg = zeros(numberNodes,3);
for i = 1:3
    currentStrain = strainExtraaa(:,:,i);
    for n = 1:numberNodes
        idx = find(n==elementNodes);
            strainAvg(n,i) = sum(currentStrain(idx))/length(currentStrain(idx));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Displacements Contours %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y]=meshgrid(0:deltaX:Lx,0:deltaY:Ly);
displacementsxx_raw=displacements(1:end/2);
displacementsyy_raw=displacements(end/2+1:end);
displacementsxx=reshape(displacementsxx_raw,numberElementsX+1,numberElementsY+1);
displacementsyy=reshape(displacementsyy_raw,numberElementsX+1,numberElementsY+1);
figure;
displacementsxx_contour = contourf(X,Y,displacementsxx','edgecolor','none');
shading interp
colorbar
title('Displacement xx Contour')
axis equal
figure;
displacementsyy_contour = contourf(X,Y,displacementsyy','edgecolor','none');
shading interp
colorbar
title('Displacement yy Contour')
axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Stress Contours %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y]=meshgrid(0:deltaX:Lx,0:deltaY:Ly);
Stressmesh1xx_raw=stressAvg(:,1);
Stressmesh1yy_raw=stressAvg(:,2);
Stressmesh1xy_raw=stressAvg(:,3);
Stressmesh1xx=reshape(Stressmesh1xx_raw,numberElementsX+1,numberElementsY+1);
Stressmesh1yy=reshape(Stressmesh1yy_raw,numberElementsX+1,numberElementsY+1);
Stressmesh1xy=reshape(Stressmesh1xy_raw,numberElementsX+1,numberElementsY+1);
figure;
Stressmesh1xx_contour = contourf(X,Y,Stressmesh1xx','edgecolor','none');
shading interp
colorbar
title('Stress σxx Contour')
axis equal
figure;
Stressmesh1yy_contour = contourf(X,Y,Stressmesh1yy','edgecolor','none');
shading interp
colorbar
title('Stress σyy Contour')
axis equal
figure;
Stressmesh1xy_contour = contourf(X,Y,Stressmesh1xy','edgecolor','none');
shading interp
colorbar
title('Stress 𝜏xy Contour')
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Strain Contours %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y]=meshgrid(0:deltaX:Lx,0:deltaY:Ly);
Strainmesh1xx_raw=strainAvg(:,1);
Strainmesh1yy_raw=strainAvg(:,2);
Strainmesh1xy_raw=strainAvg(:,3);
Strainmesh1xx=reshape(Strainmesh1xx_raw,numberElementsX+1,numberElementsY+1);
Strainmesh1yy=reshape(Strainmesh1yy_raw,numberElementsX+1,numberElementsY+1);
Strainmesh1xy=reshape(Strainmesh1xy_raw,numberElementsX+1,numberElementsY+1);
figure;
Strainmesh1xx_contour = contourf(X,Y,Strainmesh1xx','edgecolor','none');
shading interp
colorbar
title('Strain εxx Contour')
axis equal
figure;
Strainmesh1yy_contour = contourf(X,Y,Strainmesh1yy','edgecolor','none');
shading interp
colorbar
title('Strain εyy Contour')
axis equal
figure;
Strainmesh1xy_contour = contourf(X,Y,Strainmesh1xy','edgecolor','none');
shading interp
colorbar
title('Strain ɣxy Contour')
axis equal