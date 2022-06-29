function [RI_map , ARI] = construct_resistivity_index(video,path)
    %https://en.wikipedia.org/wiki/Arterial_resistivity_index
    %adapted to Resistivity application 
    disp(size(video)) ;
    figure(100)
    imagesc(video(:,:,1)) ; 
    colormap gray ; 
    M0 = mat2gray(video) ;
    M0 = reshape(M0,[size(M0,1),size(M0,2),1,size(M0,3)]) ;

    tmp = zeros(2*size(M0, 1)-1, 2*size(M0, 2)-1, size(M0,3), size(M0,4));
    for mm = 1:size(M0, 3)
        for pp = 1:size(M0, 4)
            tmp(:,:,mm,pp) = squeeze(interp2(squeeze(M0(:,:,mm,pp)), 1));
        end
    end
    M0 = tmp;
    disp('size M0') ; 
    disp(size(M0)) ; 


    %% mask std 


    path_png = fullfile(path,'Resistivity','png') ; 
    path_eps = fullfile(path,'Resistivity','eps') ; 

%     if  ~(exist(fullfile(path, sprintf('Resistivity')), 'dir'))
        mkdir(fullfile(path, sprintf('Resistivity'))) ;
        mkdir(path_png) ;
        mkdir(path_eps) ;
%     end
    

    mask = squeeze(std(M0,1,4)) ;
%     mask = normalize(mask,'range',[0 1]) ; 

%     mask / max(mask,[], 'all')   ;

    mask = imbinarize(im2gray(mask), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.4);

    figure(1)
    imagesc(mask) ; 
    colormap gray
    title('Arteries mask') ; 
    axis square
    axis off
    set(gca,'LineWidth', 2); 
    arteries = squeeze(sum(M0 .* mask, [1 2])) / sum(mask, "all");
    sze = size(M0) ; 
    filename_png = fullfile(path_png,'arteries_mask') ; 
    filename_eps = fullfile(path_eps,'arteries_mask') ; 
    print('-f1','-depsc',filename_eps) ; 
    print('-f1','-dpng',filename_png) ;
    imwrite(mat2gray(mask),strcat(filename_png,'.png'),'png') ; 
%     ARI = 1 ;  
%%
    figure(2) 
    p = plot(arteries) ; 
    p.Color = 'k' ; 
    title('Pulse wave')
    xlabel('Time (s)','FontSize',16) ; 
    xlim([0 size(video,3)])
    pbaspect([1.618 1 1]) ; 
    set(gca,'LineWidth', 2); 
    filename_png = fullfile(path_png,'pulse_wave') ; 
    filename_eps = fullfile(path_eps,'pulse_wave') ; 
    print('-f2','-depsc',filename_eps) ; 
    print('-f2','-dpng',filename_png) ;
    imwrite(mat2gray(arteries),strcat(filename_png,'.png'),'png') ; 


%%

    
% 
%     disp([sze(1), sze(2)])

    ARI = (max(arteries)-min(arteries))/max(arteries) ; 
    
    RI_map = ones(sze(1),sze(2)) ; 
    for n=1:sze(1)
        for k=1:sze(2) 
            Y = ones(sze(4),1) ;
            Y(:) = M0(n,k,1,:) ;
            RI_map(n,k) = (max(Y)-min(Y))/max(Y) ; 
        end 
    end  

    c1 = [0 0 1] ; 
    c2 = [1 1 1] ; 
    c3 = [1 0 0] ; 
    
    N = 1000 ; 
    n4 = 0.005 *N ;
    n3 = 0.045 * N;
    n2 = 0.95 * N;
    n1 = 0; 
    n0 = 0 ;  
    cmap = [linspace(0, 0, n0) ; linspace(0, 0, n0) ; linspace(1, 1,n0) ];
    cmap(:,end+1:end+n1) = [linspace(c1(1),c2(1),n1) ; linspace(c1(2),c2(2), n1) ; linspace(c1(3),c2(3),n1)] ; 
    cmap(:,end+1:end+n2) = [linspace(1, 1, n2) ; linspace(1, 1, n2) ; linspace(1, 1,n2) ];
    cmap(:,end+1:end+n3) = [linspace(c2(1), c3(1), n3) ; linspace(c2(2), c3(2), n3) ; linspace(c2(3), c3(3),n3) ]; 
    cmap(:,end+1:end+n4) = [linspace(1, 1, n4) ; linspace(0, 0, n4) ; linspace(0, 0,n4) ];
    
%     disp(size(Y)) ; 
    figure(3)
    imagesc(RI_map.*mask) ;
    colormap(cmap') ; 
    colorbar
    axis square

    title('Arterial resistivity index map')

    

    
    %%

    [Nx, Ny, ~,num_frames] = size(M0);
    M0 = squeeze(mean(M0,4)) ; 
    sat = abs(RI_map);
    tol = [0.6, 0.999];
    sat = 0.7 * imadjustn(sat, stretchlim(sat(:), tol));
    disp('size sat') ; 
    disp(size(sat)) ; 
    M0 = imadjustn(M0, stretchlim(M0(:), [0.02, 0.998]));

    img = zeros(Nx, Ny, 3, 'single');

    img1 = hsv2rgb(1 * ones(Ny, Nx), sat(:,:), M0(:,:));
%     img2 = hsv2rgb(0.66 * ones(Nx, Ny), sat(:,:), M0(:,:));
%     img(:,:,:) = img1 .* (M0_diff(:,:) > 0) + img2 .* (M0_diff(:,:) < 0);

    %FIXME : create adequate img for display with holo/show_hologram (line 260)
    img = mat2gray(img1);
    figure(150)
    imagesc(img) ; 
    title_img = strcat('Arterial resistivity index : ', num2str(ARI)) ; 
    title(title_img) ;
    axis square
    axis off 
    set(gca,'LineWidth', 2); 
    filename_png = fullfile(path_png,'AretrialResistivityIndex') ; 
    filename_eps = fullfile(path_eps,'AretrialResistivityIndex') ; 
    print('-f150','-depsc',filename_eps) ; 
%     print('-f150','-dpng',filename_png) ;
    imwrite(mat2gray(img),strcat(filename_png,'.png'),'png') ; 
    filename = fullfile(path,'Resistivity','ARI_map') ; 
    print('-f150','-dpng',filename) ; 




end