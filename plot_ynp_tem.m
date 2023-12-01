clear all, close all, clc
%% Read in USGS inverted SkyTEM results
% USGS SkyTEM Results can be downloaded here: https://www.sciencebase.gov/catalog/item/613793ced34e40dd9c0cdae3
% Bedrosian, P. et al. Airborne electromagnetic survey processed data and models data release, Yellowstone National Park, Wyoming, 2016. US Geol. Surv. ScienceBase Data Release https://doi.org/10.5066/P9LVAU7W (2021).

filein = ('C:\Users\aparseki\Dropbox\YNP\USGS_YNP_SkyTEM\AirborneElectro\YellowstoneNP_2016_Norris_ProcessedData_Models.csv'); %input path to USGS csv file here
d = importSkyTEM(filein);

%% Extract and plot locations of all data within bounds, make abbreviated database to speet up later tasks

poi = [520089 4947287]; % choose point of interest in GoogleEarth using UTM units (easting | northing)
poi_buff = 2000; %buffer distance around point of interest to extract (meters)
nde = [poi(1)+poi_buff poi(1)-poi_buff poi(2)+poi_buff poi(2)-poi_buff]; %set MAXX MINX MAXY MINY bounds within the datacube to work with, UTM meters.

sitecnt = 0; %just a counter
for ii = 1:length(d) %loop through all rows in database
    if d(ii,2)<nde(1) && d(ii,2)>nde(2) && d(ii,3)<nde(3) && d(ii,3)>nde(4) %the geometric check to see if point ii is within buffer of poi
        sitecnt = sitecnt+1; %just a counter
        Dat(sitecnt,:) = d(ii,:); %collect just the data in the node
    end
end
figure
scatter3(Dat(:,2),Dat(:,3),Dat(:,4), ones(length(Dat),1).*30,log10(Dat(:,5)),'filled'); hold on

%%  Define flightline endpoints to extract, interpolate and plot selected transect
% note: use GoogleEarth to review available flight lines and slect
% endpoints very near known flightlines (if sparse) or anywhere (if 3D)
% <=====define start/end of transect to extract on the next line=======>
SF = [518483 4947529 521568 4947534]; % start and end point of a line, UTM meters
% FORMT: ST_est, ST_nrth, FIN_est, FIN_nrth
% <=======================================================================>

P = polyfit([SF(1) SF(3)],[ SF(2) SF(4)],1); % fit a line bewtween the given coordinates
xs = linspace(min([SF(1) SF(3)]),max([SF(1) SF(3)]),100);  %populate x's for plotting at 100m interval
y = P(1).*xs+P(2); %calculate y's given x's above
scatter3(xs,y,ones(length(xs),1),ones(length(xs),1).*10,'k'); %plot the line for QAQC

%% Extract LINE data to plot from the abbreviated database

radius = 60; %how far around the node to look for other points
sitecnt = 0; % just a counter
nde = [xs' y']; %input desired line from previous section
for ii = 1:length(Dat) %loop through all rows
    for jj = 1:length(xs) %
        if Dat(ii,2)>nde(jj,1)-radius && Dat(ii,2)<nde(jj,1)+radius && Dat(ii,3)>nde(jj,2)-radius && Dat(ii,3)<nde(jj,2)+radius %check to see if point is along line of interest
            sitecnt = sitecnt+1;
            linD(sitecnt,:) = Dat(ii,:); %collect just the data in the node
        end
    end
end

%figure
%scatter3(linD(:,2),linD(:,3),linD(:,4),
%ones(length(linD),1).*50,log10(linD(:,5)),'filled'); hold on % show the line

%% Plot 2D resistivity section

figure
subplot(4,4,1:12)
wd = 20; %plotting width
for i = 1:length(linD)
    for j = 1:29 %location of all inverted RHO values
        pos = [0 linD(i,2) linD(i,4)-linD(i,j+64)  linD(i,4)-linD(i,j+94)]; % the position of each pixel (2 = eastibng, 4=ground elevation, 64 is start of tops, 94 is start of bottoms
        if linD(i,j+94)<linD(i,212)
            patch([pos(2) pos(2)+wd pos(2)+wd pos(2)],[pos(3) pos(3) pos(4) pos(4)],log10(linD(i,j+4)),'EdgeAlpha',0,'FaceAlpha',1); hold on % plot the position as a patch with color scaled to Rho (Rho starts at 5)
        else
            patch([pos(2) pos(2)+wd pos(2)+wd pos(2)],[pos(3) pos(3) pos(4) pos(4)],log10(linD(i,j+4)),'EdgeAlpha',0,'FaceAlpha',0.05); hold on % plot the position as a patch with color scaled to Rho (Rho starts at 5)
        end
    end
end

plot(linD(:,2),linD(:,4)-linD(:,212),'--k') % plot the DOI contour
xlabel('UTM easting [m]')
ylabel('elevation [m]')
caxis([.1 3.5])
colormap jet

subplot(4,4,13)
caxis([.1 3.5])
colormap jet
colorbar('location','southoutside')
axis off
text(.15,-.25,'log_1_0(resistivity)')

%scatter3(xs,y,ones(length(xs),1).*max(linD(:,4)),ones(length(xs),1).*10,'k')

set(findall(gcf,'-property','FontName'),'FontName','AvantGarde')
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 8])
print -dtiff SkyTEM_BerylSpringNorth.tiff -r300
close all
