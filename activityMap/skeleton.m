classdef skeleton
    % skeleton: Class for defining and analysing skeletons of 3D objects
    % This class allows you to take a mask of an object (eg the mushroom
    % body), manually define a skeleton by clicking on a montage of the
    % mask, and define a "starting point" (eg where stimulation was locally
    % applied). Then it will create a series of evenly spaced nodes on the
    % skeleton to divide up the object evenly into 'zones' of certain
    % distances from the starting point
    %
    % Instructions for use:
    % Draw a mask based on channel 2 of an activityMap object, 'map':
    % >> mask = activityMap.outlineObject(map.f0, 2); 
    %
    % Create the skeleton object:
    % >> skel = skeleton(mask, map.pixelCal); 
    % This will then display the mask as a montage. The user then clicks a
    % few points in order to define one branch of the skeleton. (Single
    % click to add a point; double click to finish; or press <Return> when
    % you're finished.)
    % A dialog then pops up: add another branch or done? If you want to add
    % another branch, the new branch will connect to the existing skeleton:
    % from the *first* point on the new branch to whatever point on the
    % eixsting skeleton is closest
    %
    % Define a starting point:
    % >> skel = skel.defineStartingPoint(map.f0(:,:,:,2));
    % This will display the inputted image in a montage, and you click a
    % single point to say where the local stimulation occurred (again,
    % press <Return> when finished)
    %
    % Create spaced node skeleton: user defines spacing in microns
    % >> skel = skel.createSpacedNodes(20); % 20 micron spacing
    %
    % Create Voronoi divisions according spaced node skeleton
    % >> skel = skel.createVoronoiMask();
    %
    % See the results!
    % >> skel.drawSkeleton();
    
    % WARNINGS:
    % the first 2 dimensions of pixelCal (x and y) may not be in the right
    % order. This probably doesn't matter because the microscope should be
    % calibrated to give square images in xy.
    %
    % Written by Andrew Lin and Hoger Amin
    properties
        mask % binary mask defining what voxels are in the object
        voronoiMask = [] % each voxel in the mask is assigned to one node in spacedNodes
        nodes = [] % the user-defined nodes of the skeleton
        links = [] % links between the user-defined nodes (n x 3 matrix) - 1st 2 columns are the node indices; 3rd column is the branch label
        spacedNodes = [] % nodes that are 'spacing' microns apart along the skeleton
        spacedLinks = [] % links between spaced nodes
        userStartingPoint = [] % user defines a point in space
        startingPoint = [] % the point on the skeleton closest to userStartingPoint
        startingLink = [] % which link (in links) is startingPoint on
        pixelCal % how many microns per pixel in each dimension (x y z) (1x3 matrix)
        numBranches
        splitMask = [] % each voxel in the mask is assigned to one of 2 zones
        vector = [] % vector for splitting the mask 
        dividingNode = []

        %        spacing = 20 % how far apart should the nodes in spacedNodes be

        
    end
    
    methods
        % This is the old version of entering skeleton via montage
        % This was used in Amin et al 2020. Might need to use this for
        % backward compatibility
%         function obj = skeleton(inputMask, inputPixelCal)
%             %
%             
%             obj.mask = inputMask;
%             obj.pixelCal = inputPixelCal;
%             
%             % display the mask in montage form
%             maskFigHandle= figure('Name','Click points for skeleton');
%             [montageMask, montageDims] = stackToMontage(permute(obj.mask,[2 1 4 3]));
%             imagesc(montageMask);
%             axis equal, axis tight;
%             
%             % manually enter points IN ORDER of the line
%             [x_in, y_in] = getpts(maskFigHandle);
%             
%             % Then need to turn this into x y z coordinates
%             [x,y,z] = montageXYToStackXYZ(x_in,y_in,[size(obj.mask,1) size(obj.mask,2)],montageDims);
%             % each of x, y and z is a nx1 vector where n is the number of points
%                         
%             % Now link these points together in a network:
%             pixelCoords = [x, y, z];
%             realCoords = obj.pixelsToMicrons(pixelCoords);
%             
%             branchNumber = 1;
%             obj = obj.appendGraphNodes(realCoords, branchNumber);
%             
%             % while still adding new branches
%             roiLoop = 1;
%             while(roiLoop)
%                 nextAction = menu('Choose an option','Add another branch','Done');
%                 if nextAction == 1
%                     % enter the points IN ORDER: 
%                     % the first point will be linked to the existing
%                     % skeleton!
%                     [x_in, y_in] = getpts(maskFigHandle);
%                     [x,y,z] = montageXYToStackXYZ(x_in,y_in,[size(obj.mask,1) size(obj.mask,2)],montageDims);
%                     pixelCoords = [x, y, z];
%                     realCoords = obj.pixelsToMicrons(pixelCoords);
%                     branchNumber = branchNumber + 1;
%                     obj = obj.appendGraphNodes(realCoords, branchNumber);
%                 elseif nextAction == 2
%                     roiLoop = 0;
%                 end
%             end
%             close(maskFigHandle);
%             
%         end
        
        function obj = skeleton(inputMask, inputPixelCal)
            %
            
            obj.mask = inputMask;
            obj.pixelCal = inputPixelCal;
            
            % match up voxels in the image to dimensions for display so that we can
            % display microns instead of pixels
            [x,y,z] = meshgrid(obj.pixelCal(2)*(1:size(obj.mask,2)), ...
                obj.pixelCal(1)*(1:size(obj.mask,1)), ...
                obj.pixelCal(3)*(size(obj.mask,3):-1:1));
            % generate the 3D surface outline of the volume
            figHandle = figure('Name','Mask outline and skeleton - rotate for 3d view','Color','white');
            axHandle = gca;
            isonormals(obj.mask,patch(isosurface(x,y,z,obj.mask,0),'FaceColor',[.7 .7 .8], 'EdgeColor','none'));
            % apply a light
            camlight;
            % make the surface 30% transparent
            alpha(0.3);
            axis equal;
            hold on;
            
            try
                coords = click3D(axHandle);
                realCoords = coords([2 1 3]);
            catch ME
                error('could not get the 3d point');
            end
            
            disp('Waiting for you to press n for next point, or t for terminate');

            while true
                w = waitforbuttonpress;
%                 lastCharacter = figHandl
                switch figHandle.CurrentCharacter
                    case {'n'} %'n' for 'next'
                        try
                            coords = click3D(axHandle);
                            realCoords = [realCoords; coords([2 1 3])];
                        catch ME
                            error('could not get the 3d point');
                        end
                        disp('Waiting for you to press n for next point, or t for terminate');

                    case {'t'} % 't' for 'terminate'
                        break;
                end
            end
%             realCoords = obj.pixelsToMicrons(pixelCoords);
            
            branchNumber = 1;
            obj = obj.appendGraphNodes(realCoords, branchNumber);
            
            % while still adding new branches
            % THIS CODE HAS NOT BEEN TESTED!!!!!!!!!!!!
            roiLoop = 1;
            while(roiLoop)
                nextAction = menu('Choose an option','Add another branch','Done');
                if nextAction == 1
                    % enter the points IN ORDER: 
                    % the first point will be linked to the existing
                    % skeleton!
                    try
                        coords = click3D(axHandle);
                        realCoords = coords([2 1 3]);
                    catch ME
                        error('could not get the 3d point');
                    end
                    while waitforbuttonpress
                        lastCharacter = get(figHandle,'CurrentCharacter');
                        switch lastCharacter
                            case {'n'} %'n' for 'next'
                                try
                                    coords = click3D(axHandle);
                                    realCoords = [realCoords; coords([2 1 3])];
                                catch ME
                                    error('could not get the 3d point');
                                end                                
                            case {'t'} % 't' for 'terminate'
                                break;
                        end
                    end
%                     realCoords = obj.pixelsToMicrons(pixelCoords);
                    branchNumber = branchNumber + 1;
                    obj = obj.appendGraphNodes(realCoords, branchNumber);
                elseif nextAction == 2
                    roiLoop = 0;
                end
            end
            close(figHandle);
            
        end
        
        function obj = divideMask(obj)
            % You will be asked to click 2 points. Call the first point you
            % click A, second point B.
            % The mask will be divided at point A.
            % The part of the mask in the same direction as the A->B vector
            % will be in the "first" part of the resulting split mask, that
            % is obj.splitMask(:,:,:,1) will be true for these points.
            % The other part, which is "behind" you if you face from A
            % toward B, will be in the "second" part -
            % obj.splitMask(:,:,:,2)
            %
            % match up voxels in the image to dimensions for display so that we can
            % display microns instead of pixels
            [x,y,z] = meshgrid(obj.pixelCal(2)*(1:size(obj.mask,2)), ...
                obj.pixelCal(1)*(1:size(obj.mask,1)), ...
                obj.pixelCal(3)*(size(obj.mask,3):-1:1));
            % generate the 3D surface outline of the volume
            figHandle = figure('Name','Mask outline and skeleton - rotate for 3d view','Color','white');
            axHandle = gca;
            isonormals(obj.mask,patch(isosurface(x,y,z,obj.mask,0),'FaceColor',[.7 .7 .8], 'EdgeColor','none'));
            % apply a light
            camlight;
            % make the surface 30% transparent
            alpha(0.3);
            axis equal;
            hold on;
            
            w=size(obj.mask,1);
            l=size(obj.mask,2);
            h=size(obj.mask,3);

            normalVector = click3Dline(axHandle);
            normalVector = normalVector(:, [2 1 3]); %flip x and y dimensions because of the drawing
            normalVector(:,3) = h*obj.pixelCal(3)-normalVector(:,3); %flip z dimension upside down;
            n_rel = normalVector(2,:) - normalVector(1,:); %normalVectorActualLocs(2,:) - normalVectorActualLocs(1,:); % 1x3 array
            
            % get all the points where mask is not 0
            [x_vals,y_vals,z_vals]=ind2sub([w,l,h],find(obj.mask(:)));
            % make an nx3 matrix
            points = [x_vals, y_vals, z_vals];
            numPoints = size(points,1);
            p_rel = points.*repmat(obj.pixelCal, numPoints, 1) - repmat(normalVector(1,:), numPoints, 1); %pointsActualLocs - repmat(normalVectorActualLocs(1,:), numPoints, 1); % numPoints x 3 array
            
            dotProducts = p_rel*n_rel'; 
            % each element of dotProducts will be the dot product of the
            % normal vector with the vector from the first point of the
            % normal vector to each point in the mask
%             dotProducts = zeros(size(pointsActualLocs,1),1);
%             for i=1:length(dotProducts)
%                 dotProducts(i) = dot(pointsActualLocs(i,:)-normalVectorActualLocs(1,:), normalVectorActualLocs(2,:) - normalVectorActualLocs(2,:));
%             end
            obj.vector = normalVector;
            
            obj.splitMask = false(w, l, h, 2);
            for i=1:length(dotProducts)
                if dotProducts(i)>=0
                    obj.splitMask(points(i,1), points(i,2), points(i,3), 1) = true;
                else
                    obj.splitMask(points(i,1), points(i,2), points(i,3), 2) = true;
                end
            end
                        
        end

        function obj = divideSkeleton(obj)
            % add a node in the skeleton exactly where the skeleton
            % intersects the plane defined by the normal vector saved from
            % divideMask
            
            normalVector = obj.vector;
            normalVector(:,3) = size(obj.mask,3)*obj.pixelCal(3)-normalVector(:,3); % turn the normal vector upside down
            % Note: I'm not sure why you have to turn it upside down. It
            % was turned upside down in divideMask but again I'm not sure
            % why - 10 Nov 2023

            p0 = normalVector(1,:); % the first point on the normal vector
            n_rel = normalVector(2,:) - normalVector(1,:); % the normal vector
            ds = zeros(size(obj.links,1),1);
            dLessThanLinkLength = false(size(obj.links,1),1);
            intersects = zeros(size(obj.links,1),3);
            intersectInLink = false(size(obj.links,1),1);
            for i=1:size(obj.links,1)
                % the variable names here are as on Wikipedia, https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
                % as of 9 Nov 2023
                % l0 is one point on the line
                % l is the unit vector in the direction of the line
                % d is the distance you need to travel along the line from l0 to get to the plane
                l0 = obj.nodes(obj.links(i,1)).realCoords;
                l1 = obj.nodes(obj.links(i,2)).realCoords; % l1 is the other end of the link
                l = (l1 - l0); %/norm(l1 - l0);
                % if dot(l,n_rel) == 0 % the line is parallel to the plane
                %     if dot(p0-l, n_rel)==0 % the line is in the plane
                %         d = 0.5*norm(l1-l0); % might as well say the point of intersection is halfway along the line
                %     else
                %         d = NaN; % the line never intersects the plane
                %     end
                % else
                %     d = dot(p0-l, n_rel)/dot(l,n_rel); % this equation is from Wikipedia
                % end
                % ds(i) = d;
                % dLessThanLinkLength(i) = abs(d)<norm(l1-l0);

                % read the documentation for line_plane_intersection
                [intersectPoint,rc] = line_plane_intersection(l, l0, n_rel, p0);
%                 [intersects(i,:),rc] = line_plane_intersection(l, l0, n_rel, p0);
                if rc~=1
                    rc
                    disp('Warning no unique intersection')
                    intersects(i,:) = NaN; % danger! what if we get NaNs later on?
                else
                    intersects(i,:) = intersectPoint;
                    
                    % use the triangle theorem to find if the intersecting
                    % point is on the line
                    if norm(intersects(i,:)-l0) + norm(intersects(i,:)-l1) <= norm(l1-l0)+0.001 % 0.001 allows for rounding errors
                        % norm(intersects(i,:)-l0) + norm(intersects(i,:)-l1)
                        % norm(l1-l0)+0.001
                        intersectInLink(i) = true;
                        % i
                    end
                end
            end
            % ds
            % dLessThanLinkLength
            % correctLink = find(ds>0 & dLessThanLinkLength)
            % intersects
            % intersectInLink


            % %% debugging
            % figure
            % hold on
            % % draw the nodes
            % for i=1:length(obj.nodes)
            %     x = obj.nodes(i).realCoords(1);
            %     y = obj.nodes(i).realCoords(2);
            %     z = obj.nodes(i).realCoords(3);
            % 
            %     if i==obj.dividingNode
            %         colorchoice = 'k';
            %     else
            %         colorchoice = 'r';
            %     end
            %     plot3(y,x,z,'o','Markersize',5,...
            %         'MarkerFaceColor',colorchoice,...
            %         'Color','k');
            %     text(y,x,z,strcat('..',num2str(i)));
            % 
            % end
            % for j=1:size(obj.links,1)
            %     coords(1,:) = obj.nodes(obj.links(j,1)).realCoords; % 1x3
            %     coords(2,:) = obj.nodes(obj.links(j,2)).realCoords; % 1x3
            %     % coords is now a 2x3 array
            % 
            %     cmap = colormap;
            %     color = cmap(round(obj.links(j,3)*255/obj.numBranches)+1,:);
            %     line(coords(:,2)',coords(:,1)',coords(:,3)','Color',color,'LineWidth',1);
            %     % text(mean(coords(:,2)), mean(coords(:,1)), mean(coords(:,3)),strcat('.  ',num2str(obj.links(j,3))));
            % end
            % 
            % % Use this code to draw the plane: https://uk.mathworks.com/matlabcentral/answers/1989773-how-to-draw-a-plane-perpendicular-to-a-line-and-then-generate-multiple-planes-at-regular-intervals
            % % plot3(intersection(2),intersection(1),intersection(3),'o','Markersize',5,'MarkerFaceColor','k','Color','k')
            % normalVector = obj.vector;
            % normalVector(:,3) = size(obj.mask,3)*obj.pixelCal(3)-normalVector(:,3); % turn the normal vector upside down
            % 
            % line(normalVector(:,2)',normalVector(:,1)',normalVector(:,3)','Color','k','LineWidth',4)
            % plot3(normalVector(1,2),normalVector(1,1),normalVector(1,3),'o','Markersize',5,'MarkerFaceColor','k','Color','k')
            % axis equal
            if sum(intersectInLink)~=1
                intersectInLink
                % for i=1:length(intersectInLink)
                %     if intersectInLink(i)
                %         plot3(intersects(i,2),intersects(i,1),intersects(i,3),'o','Markersize',5,'MarkerFaceColor','k','Color','k')
                %     end
                % end
                error('did not find a unique intersection between the dividing plane and the skeleton')
            end
            correctLink = find(intersectInLink);
            intersection = intersects(correctLink,:);

            % insert the new node into the skeleton
            newNodeIndex = length(obj.nodes)+1;
            obj.nodes(newNodeIndex).realCoords = intersection;
            obj.nodes(newNodeIndex).links = obj.links(correctLink,[1 2]);
            % replace the old link with the link from point 1 to the new
            % link
            node1 = obj.links(correctLink,1);
            node2 = obj.links(correctLink,2);
            linkLabel = obj.links(correctLink,3);
            obj.links(correctLink,:) = [node1, newNodeIndex, linkLabel];
            newLinkIndex = size(obj.links,1)+1;
            obj.links(newLinkIndex,:) = [newNodeIndex, node2, linkLabel];
            % link the node 1 and node 2 to the new node instead of each
            % other
            node1OrigLinks = obj.nodes(node1).links;
            node2Index = find(node1OrigLinks==node2);
            obj.nodes(node1).links(node2Index) = newNodeIndex;
            node2OrigLinks = obj.nodes(node2).links;
            node1Index = find(node2OrigLinks==node1);
            obj.nodes(node2).links(node1Index) = newNodeIndex;


            obj.dividingNode = newNodeIndex;

            % obj.drawDividedSkeleton;
            

        end

        function [] = drawDividedSkeleton(obj)
            figure
            hold on
            % draw the nodes
            for i=1:length(obj.nodes)
                x = obj.nodes(i).realCoords(1);
                y = obj.nodes(i).realCoords(2);
                z = obj.nodes(i).realCoords(3);

                if i==obj.dividingNode
                    colorchoice = 'k';
                else
                    colorchoice = 'r';
                end
                plot3(y,x,z,'o','Markersize',5,...
                    'MarkerFaceColor',colorchoice,...
                    'Color','k');
                text(y,x,z,strcat('..',num2str(i)));

            end
            for j=1:size(obj.links,1)
                coords(1,:) = obj.nodes(obj.links(j,1)).realCoords; % 1x3
                coords(2,:) = obj.nodes(obj.links(j,2)).realCoords; % 1x3
                % coords is now a 2x3 array

                cmap = colormap;
                color = cmap(round(obj.links(j,3)*63/obj.numBranches)+1,:);
                line(coords(:,2)',coords(:,1)',coords(:,3)','Color',color,'LineWidth',1);
                % text(mean(coords(:,2)), mean(coords(:,1)), mean(coords(:,3)),strcat('.  ',num2str(obj.links(j,3))));
            end
            % plot3(intersection(2),intersection(1),intersection(3),'o','Markersize',5,'MarkerFaceColor','k','Color','k')
            normalVector = obj.vector;
            normalVector(:,3) = size(obj.mask,3)*obj.pixelCal(3)-normalVector(:,3); % turn the normal vector upside down

            line(normalVector(:,2)',normalVector(:,1)',normalVector(:,3)','Color','k','LineWidth',4)
            plot3(normalVector(1,2),normalVector(1,1),normalVector(1,3),'o','Markersize',5,'MarkerFaceColor','k','Color','k')
            axis equal
        end
        
        % This is the version using a montage instead of 3d-clicking
        % This was used in Amin et al 2020. Might need to use this for
        % backward compatibility
%         function obj = defineStartingPoint(obj, image)
%             % image must be a 3d image
%             image = squeeze(image);
%             if size(image)~=size(obj.mask)
%                 error('defineStartingPoint: input image is not the same size as the mask');
%             end
%             figHandle = figure('Name', 'Click on the starting point, then press <return>');
%             [montageImage, montageDims] = stackToMontage(permute(image,[2 1 4 3]));
%             imagesc(montageImage);
%             axis equal, axis tight;
%             
%             % manually enter a point
%             [x_in, y_in] = getpts(figHandle);
%             [x, y, z] = montageXYToStackXYZ(x_in,y_in,[size(obj.mask,1) size(obj.mask,2)],montageDims);
%             
%             obj.userStartingPoint = obj.pixelsToMicrons([x, y, z]);
%             
%             close(figHandle);
%         end
%   
        function obj = defineStartingPoint(obj, varargin)
            % varargin for backward compatibility, but this code will
            % ignore the input image and just use the saved mask
            % match up voxels in the image to dimensions for display so that we can
            % display microns instead of pixels
            [x,y,z] = meshgrid(obj.pixelCal(2)*(1:size(obj.mask,2)), ...
                obj.pixelCal(1)*(1:size(obj.mask,1)), ...
                obj.pixelCal(3)*(size(obj.mask,3):-1:1));
            % generate the 3D surface outline of the volume
            figHandle = figure('Name','Mask outline and skeleton - rotate for 3d view','Color','white');
            axHandle = gca;
            isonormals(obj.mask,patch(isosurface(x,y,z,obj.mask,0),'FaceColor',[.7 .7 .8], 'EdgeColor','none'));
            % apply a light
            camlight;
            % make the surface 30% transparent
            alpha(0.3);
            axis equal;
            hold on;
            % draw the nodes
            for i=1:length(obj.nodes)
                x = obj.nodes(i).realCoords(1);
                y = obj.nodes(i).realCoords(2);
                z = obj.nodes(i).realCoords(3);
                
                plot3(y,x,z,'o','Markersize',5,...
                    'MarkerFaceColor','r',...
                    'Color','k');
                
            end
            for j=1:size(obj.links,1)
                coords(1,:) = obj.nodes(obj.links(j,1)).realCoords; % 1x3
                coords(2,:) = obj.nodes(obj.links(j,2)).realCoords; % 1x3
                % coords is now a 2x3 array
                
                line(coords(:,2)',coords(:,1)',coords(:,3)','Color','r','LineWidth',1);
                %                 text(mean(coords(:,2)), mean(coords(:,1)), mean(coords(:,3)),strcat('.  ',num2str(obj.links(j,3))));
            end
            
            try
                coords = click3D(axHandle);
                realCoords = coords([2 1 3]);
            catch ME
                error('could not get the 3d point');
            end
            
            obj.userStartingPoint = realCoords;
            
            close(figHandle);
        end
        
        
        function obj = appendGraphNodes(obj, realCoords, branchNumber)
            % link up all the inputted points in realCoords into a single linear
            % branch
            branchNodes = [];
            branchLinks = [];
            for i=1:size(realCoords,1)
                branchNodes(i).realCoords = realCoords(i,:);
                if i<2
                    prevNode = [];
                else
                    prevNode = i-1;
                end
                if i==size(realCoords,1)
                    nextNode = [];
                else
                    nextNode = i+1;
                end
                branchNodes(i).links = [prevNode nextNode];
                if ~isempty(nextNode)
                    branchLinks = [branchLinks; [i, nextNode, branchNumber]];
                end
            end
            
            % append this branch to the nodes/links of this object
            
            if ~isempty(obj.nodes) && ~isempty(obj.links)
                
                % find the point on the existing links that is closest to the FIRST node on
                % the new branch
                [branchPoint, branchLink] = obj.closestPointOnSkeleton(branchNodes(1).realCoords);
                
                % add a new node at branch point
                obj.nodes(end+1).realCoords = branchPoint;
                branchPointIndex = length(obj.nodes);
                % de-link the 2 ends of branchLink from each other and link them instead to
                % branchPoint
                node1 = obj.links(branchLink,1);
                node2 = obj.links(branchLink,2);
                oldBranchNumber = obj.links(branchLink,3);
                % delete the old link
                obj.links(branchLink,:) = [];
                % add new links to the new branch point
                obj.links = [obj.links; [node1, branchPointIndex, oldBranchNumber]; [node2, branchPointIndex, oldBranchNumber]];
                % replace the links in the nodes and link them to branchPoint
                obj.nodes(node1).links(obj.nodes(node1).links==node2) = branchPointIndex;
                obj.nodes(node2).links(obj.nodes(node2).links==node1) = branchPointIndex;
                % link new branchPoint to the old nodes
                obj.nodes(branchPointIndex).links = [node1 node2];
                
                % Now add the new branch
                % need to re-index all the node numbers in branchNodes
                for i=1:length(branchNodes)
                    branchNodes(i).links = branchNodes(i).links + length(obj.nodes);
                end
                % and re-index the node numbers in branchLinks
                branchLinks(:,[1 2]) = branchLinks(:,[1 2]) + length(obj.nodes);
                % add the links
                obj.links = [obj.links; branchLinks];
                % add the nodes
                obj.nodes = [obj.nodes, branchNodes];
                
                % link new branchPoint to the old nodes plus the first node of the branch
                obj.nodes(branchPointIndex).links = [obj.nodes(branchPointIndex).links branchPointIndex+1];
                obj.nodes(branchPointIndex+1).links = [obj.nodes(branchPointIndex+1).links branchPointIndex];
                
                % add the link between the new branch and the old branch
                obj.links = [obj.links; [branchPointIndex, (branchPointIndex+1), branchNumber]];
            else
                obj.nodes = branchNodes;
                obj.links = branchLinks;
            end
        end
        
        function obj = labelBranches(obj)
            % labels the links in the skeleton 
            % if you drew the skeleton in this order:
            % vertical lobe tip -> calyx, then add horizontal lobe on
            % then vertical lobe = 1, peduncle/calyx = 2, horizontal lobe =
            % 3 [that was the way in Amin et al. 2020]
            % in contrast if you draw the skeleton in this order:
            % calyx -> horizontal lobe tip, then add vertical lobe 
            % then peduncle/calyx = 1, horizontal lobe = 2, vertical lobe =
            % 3
            
            % find all nodes with >2 connections
            nodeLinks = {obj.nodes.links}';
            numLinksPerNode = zeros(length(nodeLinks));
            for j=1:length(nodeLinks)
                numLinksPerNode(j) = length(nodeLinks{j});
            end
            junctionIndices = find(numLinksPerNode>2);
            % initialize all branch labels
            obj.links(:,3) = 0;
            %loop through the junctions
            count = 1;
            for i=1:length(junctionIndices)
                % get the neighbors of the junction
                neighbors = nodeLinks{junctionIndices(i)};
                neighbors = sort(neighbors);
                for j=1:length(neighbors)
                    reachedEnd = false;
                    currentNode = junctionIndices(i); 
                    nextNode = neighbors(j);
                    linkIndex = obj.getLinkIndices([currentNode nextNode]);
                    if ~obj.links(linkIndex,3)
                        while ~reachedEnd
                            obj.links(linkIndex,3) = count;
                            nextNeighbors = obj.nodes(nextNode).links;
                            if length(nextNeighbors) ~= 2
                                reachedEnd = true;
                            else
                                lastNode = currentNode;
                                currentNode = nextNode;
                                nextNode = nextNeighbors(nextNeighbors~=lastNode);
                                linkIndex = obj.getLinkIndices([currentNode nextNode]);
                            end
                        end
                        count = count+1;
                    end
                end
            end

            obj.numBranches = 3;
        end

        function obj = labelBranchesSepCalyx(obj)
            % labels the links in the skeleton 
            % calyx = 1, peduncle = 2, horizontal = 3, vertical = 4
            % It assumes that you draw the skeleton in this order:
            % calyx -> horizontal lobe tip, then add vertical lobe 
            % then peduncle/calyx = 1, horizontal lobe = 2, vertical lobe =
            % 3
            
            % find all nodes with >2 connections
            nodeLinks = {obj.nodes.links}';
            numLinksPerNode = zeros(length(nodeLinks));
            for j=1:length(nodeLinks)
                numLinksPerNode(j) = length(nodeLinks{j});
            end
            junctionIndices = find(numLinksPerNode>2);
            % Treat the dividing node as a junction
            junctionIndices = [obj.dividingNode junctionIndices ];
            % initialize all branch labels
            obj.links(:,3) = 0;
            %loop through the junctions
            count = 1;
            for i=1:length(junctionIndices)
                % get the neighbors of the junction
                neighbors = nodeLinks{junctionIndices(i)};
                neighbors = sort(neighbors);
                for j=1:length(neighbors)
                    reachedEnd = false;
                    currentNode = junctionIndices(i); 
                    nextNode = neighbors(j);
                    linkIndex = obj.getLinkIndices([currentNode nextNode])
                    if ~obj.links(linkIndex,3)
                        while ~reachedEnd
                            obj.links(linkIndex,3) = count;
                            nextNeighbors = obj.nodes(nextNode).links;
                            if length(nextNeighbors) ~= 2 || ismember(nextNode, junctionIndices) %%% IS THIS GOING TO WORK??? UNKNOWN????
                                reachedEnd = true;
                            else
                                lastNode = currentNode;
                                currentNode = nextNode;
                                nextNode = nextNeighbors(nextNeighbors~=lastNode);
                                linkIndex = obj.getLinkIndices([currentNode nextNode]);
                            end
                        end
                        count = count+1;
                    end
                end
            end

            obj.numBranches = 4;
        end
        
        function result = splitBranches(obj)
            % Return an n x numPaths matrix where each column follows
            % the skeleton from the starting point, along every possible
            % path
            
        end
        
        function obj = createSpacedNodes(obj, spacing)
            obj.numBranches=3;
            disp('Creating spaced nodes');
            if isempty(obj.userStartingPoint)
                error('createSpacedNodes: you must run defineStartingPoint first');
            end
            if isscalar(spacing)
                spacing = repmat(spacing, obj.numBranches, 1);
            end
            
            [obj.startingPoint, obj.startingLink] = obj.closestPointOnSkeleton(obj.userStartingPoint);
                        
            startingNode = [];
            startingNode.realCoords = obj.startingPoint;
            startingNode.distFromStart = 0;
            startingNode.skeletonLink = obj.startingLink;
            startingNode.branchNumber = obj.links(obj.startingLink,3);
            startingNode.links = []; % links to other spaced nodes
            
            % we will construct a new array of nodes that are all evenly spaced
            % they will be labeled with the distance from the starting node
            % each node will be a struct with the fields:
            % - realCoords
            % - distance from startingNode
            % - which skeletonLink is it on
            % - which branchNumber is it on
            % - which other nodes is it connected to
            
            % the new nodes
            obj.spacedNodes = startingNode;
            obj.spacedLinks = [];
            
            skeletonNodes = obj.nodes;
            skeletonLinks = obj.links;
            
            % breadth-first search in case there are loops in the skeleton
            spacedNodeQueue = [];
            spacedNodeQueue = [spacedNodeQueue; 1];
            % the stack contains INDICES of nodes (in spacedNodes), not the nodes themselves
            
            % mark which skeletonLinks have been visited already
            % initialise with all at 0 - they will be marked as 1 as they get visited
            skeletonLinksVisited = zeros(length(skeletonLinks), 1);
            
            % key to variable names:
            % t: the starting spacedNode
            % n: the successive skeletonNodes after t
            % o: the skeletonNode after n
            distanceTolerance = 0.01;
            while ~isempty(spacedNodeQueue)
                % take the top element of spacedNodeQueue BUT DO NOT POP IT
                % OFF YET (we only remove spacedNodes when we are sure
                % you cannot create any more spacedNodes from it)
                tIndex = spacedNodeQueue(1);
                t = obj.spacedNodes(tIndex);
                % mark the link that t is on as visited
                skeletonLinksVisited(t.skeletonLink) = 1;
                
                % test if t is actually on the link that it's supposed to be on
                linkEndpoint1 = skeletonNodes(skeletonLinks(t.skeletonLink,1)).realCoords;
                linkEndpoint2 = skeletonNodes(skeletonLinks(t.skeletonLink,2)).realCoords;
                if abs(pdist2(t.realCoords,linkEndpoint1) + pdist2(t.realCoords,linkEndpoint2) ...
                        - pdist2(linkEndpoint1, linkEndpoint2)) > distanceTolerance
                    error('Node is not on the link that it is supposed to be on!');
                end
                
                createdASpacedNode = 0;
                % for both ends of the current link
                for i=1:2
                    newCoords = NaN;
                    nextSkeletonNodeIndex = skeletonLinks(t.skeletonLink,i);
                    nextSkeletonNode = skeletonNodes(nextSkeletonNodeIndex);
                    vector = nextSkeletonNode.realCoords - t.realCoords;
                    %display(strcat('norm(vector)=',num2str(norm(vector)))); % debugging delete %%%%%%%%%%%%%%%%%%%%%%%
                    unitVector = vector / norm(vector);
                    branchNumber = skeletonLinks(t.skeletonLink,3);
                    %display(strcat('branchNumber=',num2str(branchNumber))); %%%%%%%%%%%%%%%%%%%
                    if norm(vector)>=spacing(branchNumber)
                        % these two values carry down below to the
                        % commands to create a new spacedNode
                        newCoords = t.realCoords + spacing(branchNumber)*unitVector;
                        linkIndex = t.skeletonLink; % we are still on the same link
                    else
                        % DEPTH-FIRST SEARCH to find the next place to put a node at
                        % the right distance from the current node
                        % from the next node, only have remainingDistance left to
                        % search
                        remainingDistance = spacing(branchNumber) - norm(vector);
                        %disp(strcat('remainingDistance:',num2str(remainingDistance)));
                        skeletonNodeStack = [];
                        distanceStack = [];
                        lastSpacingStack = [];
                        % need to keep track of links visited on this
                        % search so that if at any point we encounter a
                        % skeletonNode with unexplored links, we will mark
                        % all links from this search as unvisited, so that
                        % we will go back to them
                        linksVisitedOnThisSearch = [];
                        % push next node and next distance onto the respective stacks
                        % only push the next node on if it's not a terminal node
                        if length(nextSkeletonNode.links)>1
                            skeletonNodeStack = [nextSkeletonNodeIndex; skeletonNodeStack];
                            distanceStack = [remainingDistance; distanceStack];
                            lastSpacingStack = [spacing(branchNumber); lastSpacingStack];
                            while ~isempty(skeletonNodeStack)
                                % nIndex <- nodeStack.pop()
                                nIndex = skeletonNodeStack(1);
                                skeletonNodeStack = skeletonNodeStack(2:end);
                                n = skeletonNodes(nIndex);
                                dist = distanceStack(1);
                                distanceStack = distanceStack(2:end);
                                lastSpacing = lastSpacingStack(1);
                                lastSpacingStack = lastSpacingStack(2:end);
                                % for all links from node n
                                skelNodesLinkedFromN = sort(n.links);
                                numLinks = length(skelNodesLinkedFromN);
                                % call getLinkIndices by horcat n (node1) and skelNodesLinkedFromN (column vector of node2)
                                linkIndices = obj.getLinkIndices([repmat(nIndex,numLinks,1) skelNodesLinkedFromN(:)]);
                                % the indices of visited match the
                                % indices of skelNodesLinkedFromN
                                % and linkIndices
                                visited = skeletonLinksVisited(linkIndices);
                                if sum(visited) < numLinks % if some of the links are unvisited
                                    % pick the first unvisited link
                                    [~,firstUnvisited] = min(visited);
                                    linkIndex = linkIndices(firstUnvisited);
                                    % mark this link as visited
                                    skeletonLinksVisited(linkIndex) = 1;
                                    % if node n still has other unvisited links after this new link has been
                                    % visited
                                    if sum(visited) < (numLinks-1)
                                        % set all previous visited links from this search to be unvisited
                                        skeletonLinksVisited(linksVisitedOnThisSearch) = 0;
                                    end
                                    % add this link to linksVisitedOnThisSearch
                                    linksVisitedOnThisSearch = [linksVisitedOnThisSearch linkIndex];
                                    oIndex = skelNodesLinkedFromN(firstUnvisited);
                                    o = skeletonNodes(oIndex);
                                    vector = o.realCoords - n.realCoords;
                                    branchNumber = obj.links(linkIndex,3);
                                    unitVector = vector / norm(vector);
                                    %disp(strcat('Now on the link between_',num2str(obj.links(linkIndex,1)),'_and_',num2str(obj.links(linkIndex,2)),'_which is branch_',num2str(branchNumber)));
                                    %disp(strcat('norm(vector)insideDFS: ',num2str(norm(vector))));
                                    %disp(strcat('dist:',num2str(dist)));
                                    thisDist = dist*spacing(branchNumber)/lastSpacing; % use thisDist in case >1 unvisited linksFromN
                                    %disp(strcat('thisDist: ',num2str(thisDist)));
                                    if norm(vector)>=thisDist %%%%%%%
                                        newCoords = n.realCoords + thisDist*unitVector;
                                        skeletonLinksVisited(linkIndex) = 1; % ensure that the link with the next node is marked visited
                                    else
                                        %display(strcat('distbeforesub ',num2str(dist)));
                                        dist = thisDist - norm(vector); %%%%%%
                                        %display(strcat('dist:', num2str(dist)));
                                        distanceStack = [dist; distanceStack];
                                        skeletonNodeStack = [oIndex; skeletonNodeStack];
                                        lastSpacingStack = [spacing(branchNumber); lastSpacingStack];
                                    end
                                %else
                                %    disp(['all links visited! from skeletonNode ', num2str(nIndex)]);
                                end
                            end
                        end
                    end
                    
                    % create the new spacedNode
                    if ~isnan(newCoords)
                        % check if you have created this one before
                        createdBefore = 0;
                        for j=1:length(obj.spacedNodes)
                            if pdist2(obj.spacedNodes(j).realCoords, newCoords) < distanceTolerance
                                createdBefore = 1;
                                break;
                            end
                        end
                        % check if you are putting a node on an existing link
                        onExistingLink = 0;
                        for j = 1:size(obj.spacedLinks,1)
                            % test if newCoords lie on any existing links
                            linkEndpoint1 = obj.spacedNodes(obj.spacedLinks(j,1)).realCoords;
                            linkEndpoint2 = obj.spacedNodes(obj.spacedLinks(j,2)).realCoords;
                            if abs(pdist2(newCoords,linkEndpoint1) + pdist2(newCoords,linkEndpoint2) ...
                                    - pdist2(linkEndpoint1, linkEndpoint2)) < distanceTolerance
                                onExistingLink = 1;
                                break;
                            end
                        end
                        if ~createdBefore && ~onExistingLink
                            newNode.realCoords = newCoords;
                            newNode.distFromStart = t.distFromStart + spacing(branchNumber);
                            newNode.skeletonLink = linkIndex;
                            newNode.links = tIndex;
                            newNode.branchNumber = skeletonLinks(linkIndex,3);
                            %display(['Created node ', num2str(tIndex), ' at distance ', num2str(newNode.distFromStart), ' on branch ', num2str(newNode.branchNumber)]);
                            obj.spacedNodes = [obj.spacedNodes; newNode];
                            % link the previous node to the new node
                            obj.spacedNodes(tIndex).links = [obj.spacedNodes(tIndex).links; length(obj.spacedNodes)];
                            obj.spacedLinks = [obj.spacedLinks; [tIndex length(obj.spacedNodes)]];
                            spacedNodeQueue = [spacedNodeQueue; length(obj.spacedNodes)];
                            createdASpacedNode = 1;
                        end
                    end

                end %for i=1:2
                if ~createdASpacedNode
                    % remove t from spacedNodeStack
                    spacedNodeQueue(spacedNodeQueue==tIndex) = [];
                end
            end % while ~isempty(spacedNodeStack)
        end % function
        
        function result = getLinkIndices(obj, nodeArray)
            % finds the indices for the links between nodes
            % nodeArray is an n x 2 array where each row is a pair of node
            % indices
            % returns an n x 1 array where the i'th element is the index of
            % the link between the nodes in the i'th row of nodeArray
            % the i'th element will be 0 if there is no such link (and will
            % print out a warning)
            % e.g.
            % nodeArray = [1 2; 3 4; 1 3];
            % getLinkIndices(nodeArray) returns the indices for the links
            % between nodes 1&2, 3&4, 1&3
            n = size(nodeArray,1);
            result = zeros(n,1);
            for i=1:n
                node1 = nodeArray(i,1);
                node2 = nodeArray(i,2);
                [~, linkIndex1] = ismember([node1 node2], obj.links(:,[1 2]), 'rows');
                [~, linkIndex2] = ismember([node2 node1], obj.links(:,[1 2]), 'rows');
                linkIndex = max(linkIndex1, linkIndex2);
                if ~linkIndex
                    disp(['getLinkIndices warning: returning 0 in element ' num2str(i)]);
                end
                result(i) = linkIndex;
            end
        end
        
        function endNodes = findEndNodes(obj)
            % Returns the indices of the nodes that have only one node
            % linked to them (i.e. "ends")
            endNodes = [];
            for i=1:length(obj.nodes)
                if length(obj.nodes(i).links)==1
                    endNodes = [endNodes i];
                end
            end
        end
        
        
        function endSpacedNodes = findEndSpacedNodes(obj)
            % Returns the indices of the spaced nodes that have only one node
            % linked to them (i.e. "ends")
            endSpacedNodes = [];
            for i=1:length(obj.spacedNodes)
                if length(obj.spacedNodes(i).links)==1
                    endSpacedNodes = [endSpacedNodes i];
                end
            end
        end
        
        function [distances, prevNodes] = getDistances(obj, startingNode)
            % This function will find the distances between the starting
            % node (the node with index number 'startingNode') and every
            % other node in the skeleton
            % Returns:
            % distances: an array whose i-th element is the distance between
            % the starting node and the node in the skeleton with index i
            % prevNodes: an array whose i-th element is the previous node
            % in the shortest path from the starting node to the node in
            % the skeleton with index i
            
            % Implements the pseudocode from the Wikipedia page about
            % Dijkstra's algorithm
            distances = Inf(length(obj.nodes), 1);
            prevNodes = zeros(length(obj.nodes), 1);
            distances(startingNode) = 0;
            Q = 1:length(obj.nodes);
            while ~isempty(Q)
                % u is the index of the node with the minimum distance 
                % but only search for the distances of nodes that are still
                % in Q
                [~, minDistIndex] = min(distances(Q));
                u = Q(minDistIndex);
                % remove u from Q
                Q = Q(Q~=u); % ie keep only values in Q that are not u
                
                % for each neighbor of u
                for i=1:length(obj.nodes(u).links)
                    % v is the index of the neighbor
                    v = obj.nodes(u).links(i);
                    % alt is the distance of node u, plus the distance
                    % between u and the neighbor
                    
                    alt = distances(u) + pdist2(obj.nodes(u).realCoords, obj.nodes(v).realCoords);
                    if (alt < distances(v)) % if we've found a shorter distance for v
                        distances(v) = alt;
                        prevNodes(v) = u;
                    end
                end
            end
        end
        
        function obj = createVoronoiMask(obj)
            anchorPoints = cell2mat({obj.spacedNodes.realCoords}');
            w=size(obj.mask,1);
            l=size(obj.mask,2);
            h=size(obj.mask,3);
            % get all the points where mask is not 0
            [x_vals,y_vals,z_vals]=ind2sub([w,l,h],find(obj.mask(:)));
            % make an nx3 matrix
            points = [x_vals, y_vals, z_vals];
            pointsActualLocs = obj.pixelsToMicrons(points);
            
            % get distances between all voxels and the 'anchor' points
            distances = pdist2(pointsActualLocs,anchorPoints);
            % find which point gives the minimum distance for each voxel
            [~,indices] = min(distances,[],2);
            
            obj.voronoiMask = false(w, l, h, size(anchorPoints,1));
            for i=1:nnz(obj.mask(:))
                obj.voronoiMask(points(i,1),points(i,2), points(i,3), indices(i)) = true;
            end
        end
        
        function obj = trimVoronoiMaskEndpoints(obj)
            % assume that for voronoiMask, the 1st and last indices of the
            % 4th dimension are the ones to trim
            % Warning: this only works for single branch skeletons!
            % NEW!! SEPT 2022 Now do it for all the spaced nodes that are
            % "end" spaced nodes
            % NOTE it does NOT work if the skeleton has only 2 nodes!!!!
            w=size(obj.voronoiMask,1);
            l=size(obj.voronoiMask,2);
            h=size(obj.voronoiMask,3);
            endNodes = obj.findEndNodes;
            for i=obj.findEndSpacedNodes %[1 size(obj.voronoiMask,4)]
                [x_vals,y_vals,z_vals]=ind2sub([w,l,h],find(reshape(obj.voronoiMask(:,:,:,i),1,[])'));
                % make an nx3 matrix
                points = [x_vals, y_vals, z_vals];
                %size(points)
                pointsActualLocs = obj.pixelsToMicrons(points);
                linkOfThisSpacedNode = obj.spacedNodes(i).skeletonLink;
                endNode = intersect(obj.links(linkOfThisSpacedNode,1:2), endNodes); % find the node on this link which is an end node
                if length(endNode)>1
                    error('this skeleton has only 2 nodes')
                end
                endPoint = obj.nodes(endNode).realCoords;
                skeletonVector = obj.spacedNodes(i).realCoords - endPoint;
%                 if i==1
%                     endPoint = obj.nodes(1).realCoords;
%                     skeletonVector = obj.nodes(2).realCoords - endPoint;
%                 elseif i==size(obj.voronoiMask,4)
%                     endPoint = obj.nodes(end).realCoords;
%                     skeletonVector = obj.nodes(end-1).realCoords - endPoint;
%                 else
%                     error('oops');
%                 end
                for j=1:size(points,1)
                     
                    vector = pointsActualLocs(j,:) - endPoint;
                    
                    if dot(vector, skeletonVector)<0
                        obj.voronoiMask(points(j,1),points(j,2),points(j,3),i) = false;
                    end
                end
            end
            
        end
        
        function [] = drawSplitSkeleton(obj, varargin)
            figure('Name','Outline and split zones - rotate for 3d view','Color','white');
            % match up voxels in the image to dimensions for display so that we can
            % display microns instead of pixels
            [x,y,z] = meshgrid(obj.pixelCal(2)*(1:size(obj.mask,2)), ...
                obj.pixelCal(1)*(1:size(obj.mask,1)), ...
                obj.pixelCal(3)*(size(obj.mask,3):-1:1));
            
            isonormals(obj.mask,patch(isosurface(x,y,z,obj.mask,0),'FaceColor',[.7 .7 .8], 'EdgeColor','none'));
            cmap = colormap;
            % if voronoiMask exists, divide up the mask with color
            % coding
            if (nargin==1)
                colors = cmap(round([1 2]*63/2)+1,:);
                for i=1:size(obj.splitMask,4)
                    isonormals(obj.splitMask(:,:,:,i),patch(isosurface(x,y,z,obj.splitMask(:,:,:,i),0),'FaceColor',colors(i,:), 'EdgeColor','none'));
                end
            else
                values = varargin{1};
                if length(values)~=size(obj.splitMask,4)
                    error('length(values) must equal # of compartments in splitMask');
                end
                maxValue = max(values(:));
                minValue = min([0 min(values(:))]);
                colors = zeros(length(values),3);
                for i=1:length(values)
                    if isnan(values(i))
                        colors(i,:)=NaN;
                    else
                        colors(i,:) = cmap(round((values(i)-minValue)*63/(maxValue-minValue))+1,:);
                    end
                end
                title(['Color coding = data, range ' num2str(minValue) ' to ' num2str(maxValue)]);
                for i=1:size(obj.splitMask,4)
                    if ~isnan(values(i))
                        isonormals(obj.splitMask(:,:,:,i),patch(isosurface(x,y,z,obj.splitMask(:,:,:,i),0),'FaceColor',colors(i,:), 'EdgeColor','none'));
                    end
                end
            end
            % apply a light
            camlight;
            % make the surface 30% transparent
            alpha(0.3);
            axis equal;
            
            % draw the vector for defining the split
            line(obj.vector(:,2)',obj.vector(:,1)',obj.vector(:,3)','Color','r','LineWidth',1);

        end
        
        function [] = drawSkeleton(obj, varargin)
            % Draw the 3d surface and skeleton
            % Usage:
            % skel.drawSkeleton() -->
            %  If spaced nodes have not been defined yet, draws the surface
            %  in gray and user-defined nodes in red
            %  If spaced nodes have been created, adds in the spaced nodes
            %  color-coded by distance from start point
            %  If the object has already been divided into Voronoi spaces
            %  by the spaced nodes, draw the surfaces not in gray but also
            %  color-coded by the distance of the relevant spaced node.
            % skel.drawSkeleton(values) -->
            %  Instead of color-coding the Voronoi spaces by distance,
            %  instead color-code them by the input parameter 'values'.
            %  length(values) must equal the number of Voronoi spaces
            %  To do this, you must have already defined the Voronoi
            %  spaces. If no Voronoi spaces have been defined, it will just
            %  draw the surface in gray.
            % skel.drawSkeleton(values,max) -->
            %  As above, but it will set the range of the color code to go
            %  from 0 to max
            %  
            figure('Name','Outline and skeleton - rotate for 3d view','Color','white');
            
            % match up voxels in the image to dimensions for display so that we can
            % display microns instead of pixels
            [x,y,z] = meshgrid(obj.pixelCal(2)*(1:size(obj.mask,2)), ...
                obj.pixelCal(1)*(1:size(obj.mask,1)), ...
                obj.pixelCal(3)*(size(obj.mask,3):-1:1));
            % generate the 3D surface outline of the volume
            if isempty(obj.voronoiMask)
                isonormals(obj.mask,patch(isosurface(x,y,z,obj.mask,0),'FaceColor',[.7 .7 .8], 'EdgeColor','none'));
                colors = zeros(length(obj.spacedNodes),3);
            else
                % comment this line out to not display the gray full mask
                isonormals(obj.mask,patch(isosurface(x,y,z,obj.mask,0),'FaceColor',[.7 .7 .8], 'EdgeColor','none'));
                cmap = colormap;
                % if voronoiMask exists, divide up the mask with color
                % coding 
                if (nargin==1)
                    values = zeros(length(obj.spacedNodes),1);
                    maxDist = max([obj.spacedNodes.distFromStart]);
                    colors = cmap(round([obj.spacedNodes.distFromStart]*(size(cmap,1)-1)/maxDist)+1,:);
                    title(['Color coding = distance, max dist = ' num2str(maxDist)]);
                    for i=1:size(obj.voronoiMask,4)
                        isonormals(obj.voronoiMask(:,:,:,i),patch(isosurface(x,y,z,obj.voronoiMask(:,:,:,i),0),'FaceColor',colors(i,:), 'EdgeColor','none'));
                    end
                else
                    values = varargin{1};
                    if length(values)~=size(obj.voronoiMask,4)
                        length(values)
                        size(obj.voronoiMask,4)
                        error('length(values) must equal # of Voronoi spaces');
                    end
                    if (nargin==3)
                        maxValue = varargin{2};
                        minValue = 0;

                    else
                        maxValue = max(values(:));
                        minValue = min(values(:));
                    end
                    for i=1:length(values)
                        if isnan(values(i))
                            colors(i,:)=NaN;
                        else
                            colors(i,:) = cmap(round((values(i)-minValue)*(size(cmap,1)-1)/(maxValue-minValue))+1,:);
                        end
                    end
                    title(['Color coding = data, range ' num2str(minValue) ' to ' num2str(maxValue)]);
                    for i=1:size(obj.voronoiMask,4)
                        if ~isnan(values(i))
                            isonormals(obj.voronoiMask(:,:,:,i),patch(isosurface(x,y,z,obj.voronoiMask(:,:,:,i),0),'FaceColor',colors(i,:), 'EdgeColor','none'));
                        end
                    end
                end
%                 if (nargin==0)
%                     maxDist = max([obj.spacedNodes.distFromStart]);
%                     for i=1:size(obj.voronoiMask,4)
%                         thisColor = cmap(round(obj.spacedNodes(i).distFromStart*63/maxDist)+1,:);
%                         isonormals(obj.voronoiMask(:,:,:,i),patch(isosurface(x,y,z,obj.voronoiMask(:,:,:,i),0),'FaceColor',thisColor, 'EdgeColor','none'));
%                     end
%                 else
%                     values = varargin{1};
%                     if length(values)~=size(obj.voronoiMask,4)
%                         error('length(values) must equal # of Voronoi spaces');
%                     end
%                     maxValue = max(values(:));
%                     for i=1:size(obj.voronoiMask,4)
%                         thisColor = cmap(round(values(i)*63/maxValue)+1,:);
%                         isonormals(obj.voronoiMask(:,:,:,i),patch(isosurface(x,y,z,obj.voronoiMask(:,:,:,i),0),'FaceColor',thisColor, 'EdgeColor','none'));
%                     end
%                 end
            end
            
            % apply a light
            camlight;
            % make the surface 30% transparent
            alpha(0.3);
            axis equal;
            
            hold on
            % draw the user-defined nodes
            for i=1:length(obj.nodes)
                x = obj.nodes(i).realCoords(1);
                y = obj.nodes(i).realCoords(2);
                z = obj.nodes(i).realCoords(3);
                
                plot3(y,x,z,'o','Markersize',5,...
                    'MarkerFaceColor','r',...
                    'Color','k');
                % text(y,x,z,strcat('..',num2str(i)));
                
            end

            % draw the skeleton
            for j=1:size(obj.links,1)
                coords(1,:) = obj.nodes(obj.links(j,1)).realCoords; % 1x3
                coords(2,:) = obj.nodes(obj.links(j,2)).realCoords; % 1x3
                % coords is now a 2x3 array
                
                line(coords(:,2)',coords(:,1)',coords(:,3)','Color','r','LineWidth',1);
                % text(mean(coords(:,2)), mean(coords(:,1)), mean(coords(:,3)),strcat('.  ',num2str(obj.links(j,3))));
            end
            
            if ~isempty(obj.spacedNodes)
                % draw the spaced nodes
                for i=1:length(obj.spacedNodes)
                    x = obj.spacedNodes(i).realCoords(1);
                    y = obj.spacedNodes(i).realCoords(2);
                    z = obj.spacedNodes(i).realCoords(3);
                    
                    % draw the node
                    if (nargin==1)
                        plot3(y,x,z,'o','Markersize',9,...
                            'MarkerFaceColor',colors(i,:),...
                            'Color','k');
                            % text(y,x,z,strcat('..',num2str(i)));
                    else
                        if ~isnan(values(i))
                            plot3(y,x,z,'o','Markersize',9,...
                                'MarkerFaceColor',colors(i,:),...
                                'Color','k');
                            %text(y,x,z,strcat(num2str(obj.spacedNodes(i).distFromStart)));
                            %text(y,x,z,strcat('..',num2str(i)));
                        end
                    end
                    
                end
            end

        end
        
        function result = pixelsToMicrons(obj, pixelCoords)
            maxZ = size(obj.mask,3);
            % turn the z dimension upside down
            pixelCoords(:,3) = maxZ - pixelCoords(:,3) + 1;
            result = pixelCoords.*repmat(obj.pixelCal,size(pixelCoords,1),1);
        end
        
        function [closestPoint, closestLink] = closestPointOnSkeleton(obj, P)
            % Find the point and link on the skeleton closest to P
            % usage:
            % result = obj.closestPointOnSkeleton(P)
            % P is an [x y z] triple in real space (not pixelspace)
            % returns an [x y z] triple (1x3) of th
            
            % start from the first node
            minDist = pdist2(obj.nodes(1).realCoords, P);
            closestPoint = obj.nodes(1).realCoords;
            closestLink = 1;
            % loop through all links
            for k = 1:size(obj.links,1)
                % the point on one end of the link
                A = obj.nodes(obj.links(k,1)).realCoords;
                % N is the vector pointing from A to the other end of the
                % link
                N = obj.nodes(obj.links(k,2)).realCoords - A;
                normN = norm(N); % the length of N
                % make the unit vector in the direction of N
                unitN = N/normN;
                % see equation at https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
                % The equation of a line can be given in vector form: X = A + t*unitN
                % in our case t is given by dot(P-A,unitN). But we don't want t to be
                % negative or to be bigger than normN or else it won't be within the line segment
                closestPointOnLine = A + min(normN,max(0,dot(P-A,unitN)))*unitN;
                distance = pdist2(closestPointOnLine, P);
                if distance < minDist
                    closestPoint = closestPointOnLine;
                    minDist = distance;
                    closestLink = k;
                end
            end
        end
        
        function result = closestSpacedNodeToUserClick(obj)
            % Asks user to click a point or points on the montage view of
            % the skeleton's mask, then returns the index or indices of the
            % spaced node(s) closest to the user-selected point(s)
            %
            % Usage:
            % result = skel.closestSpacedNodeToUserClick();
            % 
            % result will contain the indices of the closest spaced nodes
            % to the user-selected points, in the order that the user
            % clicked the points.
            %
            
            figHandle = figure('Name', 'Click on a point(s) to find the closest spaced node(s), then press <return>');
            [montageImage, montageDims] = stackToMontage(permute(obj.mask,[2 1 4 3]));
            imagesc(montageImage);
            axis equal, axis tight;
            
            % manually enter a point
            [x_in, y_in] = getpts(figHandle);
            [x, y, z] = montageXYToStackXYZ(x_in,y_in,[size(obj.mask,1) size(obj.mask,2)],montageDims);
            points = [x, y, z];
            pointsActualLocs = obj.pixelsToMicrons(points);
            
            % What is the closest spaced node to the user-entered point?
            anchorPoints = cell2mat({obj.spacedNodes.realCoords}');
            distances = pdist2(pointsActualLocs,anchorPoints);
            % find which point gives the minimum distance for each voxel
            [~,indices] = min(distances,[],2);
            result = indices;
        end
    end
        
    methods(Static)
        function [] = saveRotatingMovie(azimuthRange, azimuthStepSize, elevation, fileName)
            % saveRotatingMovie(azimuthRange, azimuthStepSize, elevation, fileName)
            % takes a figure with a 3d object and saves an .avi movie where
            % the object is viewed from a certain elevation and is rotated
            % through a range of azimuth angles, forward and backward.
            %
            % Usage:
            % map.skel.saveRotatingMovie(azimuthRange, azimuthStepSize, elevation, fileName);
            % azimuthRange: a 1x2 matrix with the endpoints of the angles
            %  over which you want to rotate, eg [-60 -120]
            % azimuthStepSize: how many degrees do you want to rotate the
            %  object for each frame of the movie?
            % elevation: the elevation angle you want to view the object
            %  from
            % fileName: a string with the file name you want to use to save
            %  the movie
            %
            % example usage:
            % map.skel.saveRotatingMovie([-120 -60], 1, 30,'rotatingMovie.avi')
            %
            % Usage notes:
            % * YOU MUST HAVE DRAWN THE SKELETON ALREADY AND HAVE THE FIGURE
            %   WITH THE SKELETON AS THE TOP WINDOW!!
            % * To decide on the right azimuthRange and elevation, rotate the
            %   figure manually and note the values displayed in the lower
            %   left corner, e.g. "Az: -120 El: 30" - note these values
            %   only appear when you are actively rotating the object.
            
            axis tight manual % this ensures that getframe() returns a consistent size
            if sign(azimuthRange(2)-azimuthRange(1)) ~= sign(azimuthStepSize)
                error('azimuthRange(2)-azimuthRange(1) must have the same sign as azimuthStepSize');
            end
            
            % forward rotation
            angles = azimuthRange(1):azimuthStepSize:azimuthRange(2);
            
            % reverse rotation
            angles = [angles azimuthRange(2):(-azimuthStepSize):azimuthRange(1)];
            
            for n = 1:length(angles)
                view(angles(n),elevation)
                % Capture the plot as an image
                frame = getframe(gcf);
                im = frame2im(frame);
                movie(:,:,:,n) = im;
            end
            
            v = VideoWriter(fileName);
            open(v);
            writeVideo(v,movie);
            close(v);
        end
        
        function result = outlineObject(image, channel, zinterval, scaleIntensity)
            % Allow user to manually outline the object in a montage view
            % usage:
            % mask = outlineObject(image, channel, zinterval, scaleIntensity);
            % Displays the specified channel of image in montage view, allow
            % user to draw multiple ROIs to define a 3D object, return a
            % binary mask (1 for in the object, 0 for outside the object)
            % Note: image must be 4D (x, y, z, channel)
            % Example usage to draw mask based on channel 2 of 4d stack:
            % mask = activityMap.outlineObject(stack, 2, 5, 10);
            % That means the input is the variable stack, we are taking
            % channel 2, every 5 z-slices, and scale the intensity so that
            % the max pixel intensity is 10x less than the actual max pixel
            % inensity. Change scaleIntensity if the image is saturated or
            % it's too faint to see your object
            channelForDisplay = image(:,:,1:zinterval:end,channel);
            numSlicesDisplayed = size(channelForDisplay,3);
            [montageForDisplay, montageDims] = stackToMontage(permute(channelForDisplay,[2 1 4 3]));

            % ask user to select the area
            fig = figure('Name', 'Draw the 3D ROI (see console for instructions)');
            min(montageForDisplay(:))
            max(montageForDisplay(:))
            montageForDisplayMin=min(montageForDisplay(:));
            montageForDisplayMax=max(montageForDisplay(:));
            imagesc(montageForDisplay, [montageForDisplayMin montageForDisplayMax/scaleIntensity]);
            axis equal, axis tight;
            maskMontage = multiROI(fig, gca);
            result = montageToStack(maskMontage,montageDims);
            % transpose image back into how Matlab stores the original data
            result = permute(result, [2 1 3]);
            result = logical(result);
            result = interpmask(result, linspace(1,numSlicesDisplayed,numSlicesDisplayed*zinterval));
            result = result(:,:,1:size(image,3));
            close(fig);
        end
    end
    
end