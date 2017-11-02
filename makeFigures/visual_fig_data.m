% script name: "visual_fig_data"

% parameters
L = 24;
save_it = 0;
step_size = 4;
num_of_pics = 3;
sigma_values = [0, .1, 2];
color_set    = ['r','b','m'];

% data
x   = [zeros(L/2-step_size,1); ones(4,1); -ones(4,1); zeros(L/2-step_size,1)];
rho = rand(L,1);
rho = rho/sum(rho);
X_base      = generate_observations(x, num_of_pics, 0, rho );
X_base(:,1) = x;

% drawing
figure;
for j=1:num_of_pics
    
    sigma = sigma_values(j);
    X = X_base + sigma*randn(L,num_of_pics);

    for k=1:num_of_pics
        ind = j + (k-1)*num_of_pics;
        h = subplot(num_of_pics, num_of_pics, ind);
        plot(X(:,k),color_set(j),'LineWidth',3);
        if (k==3)  %,sigma>0)
            mytext = ['\sigma = ',num2str(sigma)];
            xlabel(mytext);
           % annotation('textbox', [0, 1, 0, 0], 'string', mytext)
           % text(-12.5,0,mytext,'FontSize',20)
        end
        set(gca,'xtick',[])
        %if k>1
         %   set(gca,'ytick',[]);
        %end
        yl = max(max(abs(X(:))),1.5);
        ylim([-yl,yl])
        xlim([1,L]);
        set(gca,'FontSize',20)   
        p = get(h, 'pos');
       % p(3) = p(3) + 0.01;  %width, 3 precent more
        p(4) = p(4) + 0.03;  %height, 3 precent more
        set(h, 'pos', p);
    end
end


if save_it
    folder_name = 'VisualFigure1';  
    mkdir(folder_name);
    cd(folder_name);
    name_it = 'data_visualization';
    saveas(gcf, name_it, 'fig');
    print('-depsc2', name_it);
    cd '../'
end



