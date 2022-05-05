clear all, close all,clc

%% Calendering pressure 
% AL 28/04/22

data_dir = "C:\Users\Axel\Documents\DEM\results\electrode_calendaring\bt065N500";

surface_force_data = readlines(data_dir+"\surface_forces.dou","EmptyLineRule","skip","WhiteSpaceRule","trim");
surface_force_data = split(surface_force_data,', ');

surface_position_data = readlines(data_dir+"\surface_positions.dou","WhiteSpaceRule","trim","EmptyLineRule","skip");
surface_position_data = split(surface_position_data,', ');
id_index_surface_position = [];
for i = 1:length(surface_position_data(1,:))
    surface_position_data(1,i);
    if startsWith(surface_position_data(1,i),'ID') 
        surface_position_data(1,i);
        id_index_surface_position = [id_index_surface_position i]; 
    end
end

surface_types = extractAfter(surface_position_data(1,id_index_surface_position+1),5);

surface_indices = [];
for i = 1:length(surface_types)
        if surface_types(i) == "PointSurface"
           surface_indices = [surface_indices i]; 
        end
end

bottom_surface_force = str2double(surface_force_data(:,surface_indices(id_index_surface_position(1))+1));

periodic_BC_data = readlines(data_dir+"/periodic_bc.dou","EmptyLineRule","skip","WhiteSpaceRule","trim");
periodic_BC_data = split(periodic_BC_data,', ');
periodic_BC_x_min = str2double(periodic_BC_data(:,2));
periodic_BC_x_max = str2double(periodic_BC_data(:,3));
periodic_BC_y_min = str2double(periodic_BC_data(:,4));
periodic_BC_y_max = str2double(periodic_BC_data(:,5));
x_side_leght = periodic_BC_x_max-periodic_BC_x_min;
y_side_leght = periodic_BC_y_max-periodic_BC_y_min;

calendering_pressure_ = bottom_surface_force./(x_side_leght.*y_side_leght);

simulation_time = str2double(periodic_BC_data(:,1));
calendaring_surface_poition = str2double(surface_position_data(:,end-1));

figure
plot(simulation_time,calendering_pressure_,'-','LineWidth',2);
hold on
xlabel('Simulation time (s)')
ylabel('Calendering pressure (Pa)')
%xlim([0,3.35])
yyaxis right
ylabel('Calendering surface height (m)')
plot(simulation_time,calendaring_surface_poition,'-','LineWidth',1.2);
legend('Calendering pressure','Calendering surface height','location','Northeast')
grid on 
set(gca,'FontWeight','bold','FontSize',13)

figure
plot(calendaring_surface_poition,calendering_pressure_,'-','LineWidth',1.2);
%xlim([1.11,2])
xlabel('Calendering surface height (m)')
ylabel('Calendering pressure (Pa)')
grid on
set(gca,'FontWeight','bold','FontSize',13)


