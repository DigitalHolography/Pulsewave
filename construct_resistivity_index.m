function [RI_map , ARI] = construct_resistivity_index(video)
    %https://en.wikipedia.org/wiki/Arterial_resistivity_index
    %adapted to Resistivity application 
    disp(size(video)) ;
    figure(100)
    imagesc(video(:,:,1)) ; 
    colormap gray ; 
    M0 = mat2gray(video) ;
    M0 = reshape(M0,[size(M0,1),size(M0,2),1,size(M0,3)]) ;
    figure(10)
    imagesc(M0(:,:,1)) ; 
    tmp = zeros(2*size(M0, 1)-1, 2*size(M0, 2)-1, size(M0,3), size(M0,4));
    for mm = 1:size(M0, 3)
        for pp = 1:size(M0, 4)
            tmp(:,:,mm,pp) = squeeze(interp2(squeeze(M0(:,:,mm,pp)), 1));
        end
    end
    M0 = tmp;
    M0(1:10,1,1,1) ;
%     disp('size M0') ; 
%     M0(1:3,1:3,1,1:3) 
    mask = squeeze(std(M0,1,4)) ;
%     mask = normalize(mask,'range',[0 1]) ; 

%     mask / max(mask,[], 'all')   ;

    mask = imbinarize(im2gray(mask), 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.4);

    figure(1)
    imagesc(mask) ; 
    title('mask') ; 
    arteries = squeeze(sum(M0 .* mask, [1 2])) / sum(mask, "all");
    sze = size(M0) ; 
    ARI = 1 ;  

    figure(2) 
    plot(arteries) ; 
    title('Arteries flow')
    ARI = (max(arteries)-min(arteries))/max(arteries) ; 
% 
%     disp([sze(1), sze(2)])
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
    
    n4 = 65 ;
    n3 = 50 ; 
    n2 = 0 ;
    n1 = 120 ; 
    n0 = 75;  
    cmap = [linspace(0, 0, n0) ; linspace(0, 0, n0) ; linspace(1, 1,n0) ];
    cmap(:,end+1:end+n1) = [linspace(c1(1),c2(1),n1) ; linspace(c1(2),c2(2), n1) ; linspace(c1(3),c2(3),n1)] ; 
    cmap(:,end+1:end+n2) = [linspace(1, 1, n2) ; linspace(1, 1, n2) ; linspace(1, 1,n2) ];
    cmap(:,end+1:end+n3) = [linspace(c2(1), c3(1), n3) ; linspace(c2(2), c3(2), n3) ; linspace(c2(3), c3(3),n3) ]; 
    cmap(:,end+1:end+n4) = [linspace(1, 1, n4) ; linspace(0, 0, n4) ; linspace(0, 0,n4) ];
    
%     disp(size(Y)) ; 
    figure(3)
    imagesc(RI_map) ;
    colormap(cmap') ; 
    axis square
    title('Resistivity index map')
    
% %     my_map = [1 1 1 ; 0.66 0 1 ; 0.33 0 1; 0 0 1 ; 0 0 0.66 ; 0 0 0.33 ; 0.33 0 0 ; 0.66 0 0 ; 1 0 0 ;  1 0 0.33 ; 1 0 0.66 ] ;
% %     colormap(my_map) ; 
%     detrend_M0 = detrend(M0) ; 
%     figure(1) 
%     avg_M0 = squeeze(mean(M0, 4));
%     avg_M0 = mat2gray(avg_M0);
% 
%     
%     std_M0 = squeeze(std(M0, 1, 4)); 
%     std_M0 = mat2gray(std_M0);
% 
% 
%     
%     max_uint16 = 65535; 
%     v_diast = avg_M0 - (std_M0 / 2); 
%     v_syst = avg_M0 + (std_M0 / 2); 
%     RI = uint16(max_uint16 * (1 - ((v_syst - v_diast) ./ v_syst)));
%     
% %     RI = ind2rgb(RI, colormap_redblue(65536));
% % 
% %     figure(1)
% %     R = flip(RI,1) ;  
% % 
% %     imagesc(R) ;
% 
% %     colormap default
% 
% %     figure(2)
% %     RI(:,:,1) = RI(:,:,1)-10000 ; 
% %     imagesc(RI) ; 
% 
% 
%     levels = zeros(67,3) ; 
%     white_levels = ones(5,3) ; 
% %     for n = 1:10 
% %         white_levels(n,:) = 0.5+n*0.05 ; 
% %     end 
%     levels(1:5,:) = white_levels ; 
%     w2r_levels = ones(21,3) ; 
%     for n = 1:21 
%         w2r_levels(n,2) = 1-((n-1)*0.05) ; 
%         w2r_levels(n,3) = 1-((n-1)*0.05) ; 
%     end 
%     levels(6:26,:) = double(w2r_levels) ; 
% %     red_levels = zeros(10,3) ; 
% %     for n = 1:10
% %         red_levels(n,1) = 1-n*0.1 ;
% %     end 
% %     levels(12:21,:) = red_levels ; 
%     r2b_levels = zeros(20,3) ;
%     for n = 1:20
%         r2b_levels(n,1) = 1-n*0.05 ;
%         r2b_levels(n,3) = n*0.05 ;
%     end 
%     levels(27:46,:) = r2b_levels ; 
%     b2r_levels = zeros(21,3) ; 
%     for n = 1:21 
%         b2r_levels(n,3) = 1-((n-1)*0.05) ; 
%     end 
%     levels(47:67,:) = b2r_levels ; 
% %     blue_levels = zeros(10,3) ; 
% %     for n = 1:10
% %         blue_levels(n,3) = 1-n*0.1 ;
% %     end 
% %     levels(32:41,:) = double(blue_levels) ;
%     disp(levels) ; 
%     
% figure(2)
% %     RI(:,:,1) = RI(:,:,1)-10000 ; 
%     imagesc(RI) ; 
%     colormap(levels) ;
%     axis square
% %     colorbar
%     RI = ind2rgb(RI, colormap_redblue(65536));
%     my_map = [1 1 1 ; 0.66 0 1 ; 0.33 0 1; 0 0 1 ; 0 0 0.66 ; 0 0 0.33 ; 0.33 0 0 ; 0.66 0 0 ; 1 0 0 ;  1 0 0.33 ; 1 0 0.66 ] ;
%     figure(3) 
%     RI = flip(RI,1) ; 
%     imagesc(RI) ; 
% 


end