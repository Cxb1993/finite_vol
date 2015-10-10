%--------------------------------------------------------------------------
% MAE250: Final Project
% Written by: David Joseph Vasko
%
% Oblique Shock Reflection
%-------------------------------------------------------------------------
%Start fresh by closing and clearing everything
    close all; 
    clear all; 
    clc;
%bools  used to determin which figures we want
    generalsolver = true;
        saveher = false;
    timecomplexity = false;
    grids = false;
    zero_dis_plots = false;
    optimal_dis = false;
    effect_of_v = false;
    exact_solution = false;
    converger = true;
    
%------------------General solver------------------------------------------
    if generalsolver
        %SET GRID SCALE and DISSAPATION
            GS = 4;
            dis = 0.85;%0.35 for original method 0.85 for dis4
        %Generate our grid ad plot it;
            [x,y,X,Y] = gen_grid(GS);
             plot_grid(x,y,X,Y);
        %Generate the exact solution
            [Cexact,shocks] = exact_sol(x,y,X,Y);
        %MacCormack
            [C,fs,DELTA,STEPS] = macFAST(x,y,X,Y,dis,shocks);
     if true
        %Change Figure Names and Save
%         figure(fs(1)); 
%         title(['Mach Number using MacCormacks Method: Dissipation'...
%             ' = 0.85']);
%         saveas(fs(1),['C:\Users\David\Documents\Latex\250Project\optdis\'...
%             'mach8x' '.png'])
%         figure(fs(2)); 
%         title(['Pressure using MacCormacks Method: Dissipation'...
%             ' = 0.85']);
%         saveas(fs(2),['C:\Users\David\Documents\Latex\250Project\optdis\'...
%             'pres8x' '.png'])
        figure(fs(3)); 
        title(['Density using MacCormacks Method: Dissipation'...
            ' = 0.85']);
        saveas(fs(3),['C:\Users\David\Documents\Latex\250Project\optdis\'...
            'dens1x' '.png'])
        %Compare to exact
     end
            [~,fs] = compareToExact(Cexact,C,x);
     if saveher
        %Change Figure Names and Save
        figure(fs(1)); 
            title(['Exact and Numerical Mach Number for the Top'...
                ' Row: Dissapation = 0.85']);
        saveas(fs(1),['C:\Users\David\Documents\Latex\250Project\optdis\'...
            'topmach4' '.png'])
        figure(fs(2)); 
            title(['Exact and Numerical Mach Number for the Middle'...
                ' Row: Dissapation = 0.85']);
        saveas(fs(2),['C:\Users\David\Documents\Latex\250Project\optdis\'...
            'midmach4' '.png'])
        figure(fs(3)); 
            title(['Exact and Numerical Mach Number for the Bottom'...
                ' Row: Dissapation = 0.85']);
        saveas(fs(3),['C:\Users\David\Documents\Latex\250Project\optdis\'...
            'botmach4' '.png'])
     end
    end
%--------------------------------------------------------------------------
%-------------------------EXACT SOLUTION-----------------------------------
if exact_solution
        GS = 4;
   %Generate our grid ad plot it;
        [x,y,X,Y] = gen_grid(GS);   
    %Generate the exact solution
        [Cexact,shocks,fs] = exact_sol(x,y,X,Y);
    %Change Figure Names and Save
        figure(fs(1)); 
        title('Exact Solution for Mach Number: 4x Grid');
        saveas(fs(1),['C:\Users\David\Documents\Latex\250Project\exact\'...
            'mach4x' '.png'])
        figure(fs(2)); 
        title('Exact Solution for Pressure: 4x Grid');
        saveas(fs(2),['C:\Users\David\Documents\Latex\250Project\exact\'...
            'pres4x' '.png'])
        figure(fs(3)); 
        title('Exact Solution for Density: 4x Grid');
        saveas(fs(3),['C:\Users\David\Documents\Latex\250Project\exact\'...
            'dens4x' '.png'])
end
%------------------Zero Dissapation Plots----------------------------------
    if zero_dis_plots
        %Generate our grid ad plot it;
            [x,y,X,Y] = gen_grid(1);
        %Generate the exact solution
            [Cexact,shocks] = exact_sol(x,y,X,Y,false);
        %MacCormack
        [C,fs] = macFAST(x,y,X,Y,0,shocks,true,1,500);
%             figure(fs(1)); 
%             title(['Mach Number using MacCormacks Method: Dissipation'...
%                 ' = 0']);
%             saveas(fs(1),['C:\Users\David\Documents\Latex\250Project\'...
%                 'zero\mach1x' '.png'])
%             figure(fs(2)); 
%             title(['Pressure using MacCormacks Method: Dissipation'...
%                 ' = 0']);
%             saveas(fs(2),['C:\Users\David\Documents\Latex\250Project\'...
%                 'zero\pres1x' '.png']) 
%             figure(fs(3)); 
%             title(['Density using MacCormacks Method: Dissipation'...
%                 ' = 0']);
%         saveas(fs(3),['C:\Users\David\Documents\Latex\250Project\'...
%             'zero\dens1x' '.png']) 
        [~,fs]=compareToExact(Cexact,C,x);
            figure(fs(2)); 
            title(['Exact and Numerical Density for the Middle'...
                ' Row: Dissapation = 0']);
            saveas(fs(2),['C:\Users\David\Documents\Latex\250Project\'...
                'zero\middens1' '.png'])
            figure(fs(1)); 
            title(['Exact and Numerical Density for the Top'...
                ' Row: Dissapation = 0']);
            saveas(fs(1),['C:\Users\David\Documents\Latex\250Project\'...
                'zero\topdens1' '.png']) 
            figure(fs(3)); 
            title(['Exact and Numerical Density for the Bottom'...
                ' Row: Dissapation = 0']);
            saveas(fs(3),['C:\Users\David\Documents\Latex\250Project\'...
                'zero\botdens1' '.png']) 
    end
%--------------------------OPTIMAL DISSAPATION-----------------------------
    if optimal_dis
       %SET GRID SCALE and DISSAPATION
            GS = 2;
       %Generate our grid ad plot it;
            [x,y,X,Y] = gen_grid(GS);
       %Generate the exact solution
            [Cexact,shocks] = exact_sol(x,y,X,Y,false);
       %loop lists
            DIS = 0.1:0.05:3;
            bulk_errors = zeros(size(DIS));
       for i=1:length(DIS)
            dis = DIS(i);
       %MacCormack
            C = macFAST(x,y,X,Y,dis,shocks,false,1);
            bulk_errors(i) = compareToExact(Cexact,C,x,false); 
       end
       %plot the data
       f = figure; hold on;
       title('Bulk Error as a Function of Dissapation');
       xlabel('Artificial Dissapation'); ylabel('Bulk Error');
       plot(DIS,bulk_errors);
       saveas(f,['C:\Users\David\Documents\Latex\250Project\'...
            'EvsD6' '.png'])
    end
%--------------------------EFFET OF COURANT-----------------------------
    if effect_of_v
       %SET GRID SCALE and DISSAPATION
            GS = 2;
       %Generate our grid ad plot it;
            [x,y,X,Y] = gen_grid(GS);
       %Generate the exact solution
            [Cexact,shocks] = exact_sol(x,y,X,Y,false);
       %loop lists
           COURANT = 0.1:0.025:1;
            bulk_errors = zeros(size(COURANT));
       for i=1:length(COURANT)
            courant = COURANT(i);
       %MacCormack
            C = macFAST(x,y,X,Y,0.85,shocks,false,courant);
            bulk_errors(i) = compareToExact(Cexact,C,x,false); 
       end
       %plot the data
       f = figure; hold on;
       title('Bulk Error as a Function of Courant Number');
       xlabel('Courant Number'); ylabel('Bulk Error');
       plot(COURANT,bulk_errors);
       saveas(f,['C:\Users\David\Documents\Latex\250Project\'...
            'vvsD2' '.png'])
    end


%-----------------------------GRIDS----------------------------------------
    if grids
            [x,y,X,Y] = gen_grid(.4);
        f = plot_grid(x,y,X,Y);
        saveas(f,['C:\Users\David\Documents\Latex\250Project\'...
            'shape' '.png'])
            [x,y,X,Y] = gen_grid(1);
        f = plot_grid(x,y,X,Y);
        saveas(f,['C:\Users\David\Documents\Latex\250Project\'...
            'grid1' '.png'])
            [x,y,X,Y] = gen_grid(2);
        f = plot_grid(x,y,X,Y,false);
        saveas(f,['C:\Users\David\Documents\Latex\250Project\'...
            'grid2' '.png'])
            [x,y,X,Y] = gen_grid(4);
        f = plot_grid(x,y,X,Y,false);
        saveas(f,['C:\Users\David\Documents\Latex\250Project\'...
            'grid4' '.png'])
    end
%-------------------TIME COMPLEXITY----------------------------------------
    %Find Time complexity and plot 
    if timecomplexity
        T = 50:25:2000;
        times = zeros(2,length(T));
        for i=1:length(T)
              TAU = T(i);
              disp('TAU');
              if TAU <= 500
                  tic();
                  macdis(x,y,X,Y,1,shocks,TAU);
                  times(1,i) = toc();
              end
              tic();
              macFAST(x,y,X,Y,1,shocks,false,1,TAU);
              times(2,i) = toc();
        end
        T = T';times = times';
    %Plot time complexity
        f = figure; hold on;
        plot(T,times(:,2),'.b','MarkerSize',15); 
        plot(T(1:19),times(1:19,1),'.r','MarkerSize',15);
        title('Time Complexity: Time Steps');
        xlabel('Time Steps'); ylabel('Time (s)');
        legend('Optimized Solution', 'Baseline');   
        saveas(f,['C:\Users\David\Documents\Latex\250Project\'...
            'TCTS' '.png'])
    %Find Time complexity of grid size - points obtained by manual changing of
    %paramaters
        gs = [.2 .4 .6 .8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 4 6 8 10 12];
        Ot = [2.45 2.62 2.93 3.18 3.46 3.84 4.25 5.20 5.80 6.29...
            8.54 18.88 37.5724 79.70 142.26 216.57];
        Bt = [2.89 9.72 21.46 34.43 50.13 72.76 107.07 130.21...
            162.81 194.12 267.90];
        f = figure; hold on;
        plot(gs,Ot,'.b','MarkerSize', 15); 
        plot(gs(1:11),Bt,'.r','MarkerSize',15);
        title('Time Complexity: Grid Scale');
        xlabel('Scale Factor'); ylabel('Time (s)');
        legend('Optimized Solution', 'Baseline');
        saveas(f,['C:\Users\David\Documents\Latex\250Project\'...
            'TCGS' '.png'])
    end
             

