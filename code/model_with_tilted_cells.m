clear
close all
set(0,'DefaultAxesFontSize',18,'DefaultAxesFontName','CMR10')

% Specify the path to save data
path_to_save = '../model_outputs'
if ~exist(path_to_save, 'dir')
	mkdir(path_to_save);
end
%% Parameters to change

Lx_vector =[2 3 4];
Lx_name_vector = {'2', '3', '4'};

n_simulations = 10; % how many simulations per identical case

Lx_limit = 5; Ly_limit = 3; % simulation box limits
Ly = 21; 

N_iterations = 1000000; % how many iterations each simulation runs before it stops

% LJ parameter 
alpha = 6; 
kT_vector = logspace(4,-2,N_iterations);

%% Parameters - no need to change
d = 1; dy = d*sin(pi/3)*2; % Size of system

% list of all the tilt angles
angle_rotation_vector = [ 3.75 7.5 11.25] * pi / 180;  
angle_name_vector = {'7o5','15','22o5'};

% number of points inside the gap for each case
N_inside_density = [ 49 49 49
    73 72 66
    96 97 96 ];

initial_randomness = d/10;

variation_N = [1 1 1 1 ];
size_variation = [1 1 1 1];
pot_temp = zeros([],2);

%% Main loop 

for index_angle = 1: length(angle_rotation_vector) %%%%% LOOP IN ANGLE
    % Get the angle for this one
    angle_rotation = angle_rotation_vector( index_angle );
    angle_name = angle_name_vector{ index_angle };

    % Calculate how many particles for that angle and each Lx
    N_inside_target_now = round( N_inside_density( : ,  index_angle ));

    for i=1:length( Lx_vector )
        N_inside_target_matrix(i,:) = [N_inside_target_now(i)-variation_N(i): size_variation(i) :N_inside_target_now(i)+variation_N(i)/2]
    end

    for index_Lx = 1:length(Lx_vector) %%%%% LOOP IN LX

        Lx = Lx_vector(index_Lx);
        Lx_name = Lx_name_vector{index_Lx};

        for index_N_target_inside = 1:size(N_inside_target_matrix,2) %%%%%% LOOP ON HOW MANY PARTICLES WE WANT
            N_inside_target = N_inside_target_matrix(index_Lx,index_N_target_inside);

            for index_simulations = 1:n_simulations %%%%% LOOP IN HOW MANY IDENTICAL SIMULATIONS

                filename_to_save = ['N' num2str(N_iterations) '2-1-22_angle' angle_name '_Ly' num2str(Ly) '_Lx' Lx_name '_Ninside' num2str(N_inside_target) '_attempt' num2str( index_simulations ) ];

                % Prepare the initial points.
                % First get a bunch of points to rotate
                nx = 20; % On each side
                ny = round(Ly/dy); extra_y_nodes = 5;

                number_initial_nodes = (-extra_y_nodes:ny+extra_y_nodes)';
                coordinates_original = [ zeros(length(number_initial_nodes),1) -number_initial_nodes*dy  ;
                    ones(length(number_initial_nodes),1)*d/2 -(number_initial_nodes+0.5)*dy  ;];
                coordinates = [];

                for i=-nx:nx
                    coordinates = [coordinates;
                        coordinates_original(:,1)+d*i coordinates_original(:,2) ];
                end

                % First the left side
                coordinates_left = coordinates;

                % Keep only the points that we want
                index_final_left = find( ( coordinates_left(:,1) <= 0 ).* ( coordinates_left(:,1) >= -Lx_limit ) .* ...
                ( coordinates_left(:,2) >= -Ly-Ly_limit ).* ( coordinates_left(:,2) <= Ly_limit ));

                radius = sqrt( sum( coordinates.^2 , 2) ); angle = atan2( coordinates(:,2) , coordinates(:,1) );
                angle_new = angle + angle_rotation; coordinates_new = [ radius.*cos(angle_new) radius.*sin(angle_new) ];

                % Now the right side
                coordinates_right = [ radius.*cos(angle_new) radius.*sin(angle_new) ];

                % Keep only the points that we want
                index_final = find( ( coordinates_new(:,1) <= 0 ).* ( coordinates_new(:,1) >= -Lx_limit ) .* ...
                    ( coordinates_new(:,2) >= -Ly-Ly_limit ).* ( coordinates_new(:,2) <= Ly_limit ));

                % Duplicate to get the other side
                points_boundary = coordinates_new(index_final,:);
                points_boundary = [ points_boundary;
                    Lx-points_boundary(:,1) points_boundary(:,2) ];

                % Now we get the points inside
                n_x = round(Lx/d);
                n_y = ceil(N_inside_target / n_x);
                
                % Add a bit of randomness
                points_x = linspace(d/2 , Lx-d/2 , n_x)  + ( rand( 1 , n_x ) - 0.5 )*initial_randomness;
                points_y = linspace(-d/2 , -Ly+d/2 , n_y) + ( rand( 1 , n_y ) - 0.5 )*initial_randomness;

                points_inside = [];
                for i=1:length(points_x)
                    for j=1:length(points_y)
                        points_inside = [ points_inside; points_x(i) points_y(j) ];
                    end
                end
                
                points_inside = points_inside(1:N_inside_target,:);
                % Get the minimization running
                n_boundary = length( points_boundary);
                % points_inside = points_inside(1:N_goal,:);
                N = length(points_inside);

                x = [points_boundary(:,1); points_inside(:,1)];
                y = [points_boundary(:,2); points_inside(:,2)];

                %%%%%create a column that says whether or not the points are inside or on boundary

                iter=0;
                inside_or_outside=ones((N+n_boundary),1);
                for iter=1:n_boundary
                    inside_or_outside(iter,1)=0;
                end

                TRI = delaunay(x,y);
                potential_initial = calculate_potential_distance_neighbor_LJ_alpha( x , y , n_boundary , N , alpha , d );

                %% Minimization
                
                n_changes = 0;
                potential_evolution = zeros( N_iterations , 1);
                potential_now = potential_initial;
                for i=1:N_iterations
                    kT = kT_vector(i);
                    index_picked = ceil(rand*N) + n_boundary;
                    delta_x = (rand-0.5)*2 ; %/displacement_factor(i);
                    delta_y = (rand-0.5)*2 ; %/displacement_factor(i);

                    x_point_now = x(index_picked) + delta_x;            y_point_now = y(index_picked) + delta_y;
                    x_now = x; x_now( index_picked ) = x_point_now;     y_now = y; y_now( index_picked ) = y_point_now;

                    if ( ( x_point_now < Lx )*( x_point_now > 0 )*( y_point_now >= -Ly )*( y_point_now <= 0 ))
                        potential_new = calculate_potential_distance_neighbor_LJ_alpha( x_now , y_now , n_boundary , N , alpha , d );

                        delta_potential = potential_new - potential_now;

                        if delta_potential < 0 || rand < exp(-delta_potential/kT)
                            x = x_now;    y = y_now;
                            potential_now = potential_new;
                            n_changes = n_changes + 1;

                        end

                    end

                    potential_evolution(i) = potential_now;
                    
                end

                filename_to_save_points = ['N' num2str(N_iterations) '2-1-22_angle' angle_name '_Ly' num2str(Ly) '_Lx' Lx_name '_Ninside' num2str(N_inside_target) '_attempt' num2str( index_simulations ) '_points' '.csv'];
                points_file = [path_to_save filename_to_save_points];
                make_final_plot_called(x,y,TRI,n_boundary,points_file,inside_or_outside);
            end
        end
    end
end

%% Function for LJ Potential Calculation

function potential = calculate_potential_distance_neighbor_LJ_alpha( x , y , n_boundary , N , alpha ,  d );

TRI = delaunay(x , y );
potential = 0;

all_pairs = union( ([TRI(:,1:2) ; TRI(:,[1 3]) ]) , TRI(:,2:3) , 'rows' );

all_pairs_inside = all_pairs( find( sum( ( all_pairs > n_boundary ) ,2)),:);
all_pairs_sorted = unique( sort( all_pairs_inside , 2 ) , 'rows' );

distances = ( x( all_pairs_sorted(:,1) ) - x(all_pairs_sorted(:,2)) ).^2 + ( y( all_pairs_sorted(:,1) ) - y(all_pairs_sorted(:,2)) ).^2;

sigma = d/(2^(1/alpha));

potential = sum(  ( sigma./ distances  ).^(2*alpha) - ( sigma./ distances ).^alpha );

end

%% Function for saving CSV files

function [allpoints] = make_final_plot_called(x,y,TRI,n_boundary,filename_to_save_cvs, inside_or_outside)

potential = 0;

all_pairs = union( ([TRI(:,1:2) ; TRI(:,[1 3]) ]) , TRI(:,2:3) , 'rows' );

all_pairs_inside = all_pairs( find( sum( ( all_pairs > n_boundary ) ,2)),:);
all_pairs_sorted = unique( sort( all_pairs_inside , 2 ) , 'rows' );

allpoints=cat(2,x,y,inside_or_outside);
length(allpoints);
csvwrite(filename_to_save_cvs,allpoints)

end
