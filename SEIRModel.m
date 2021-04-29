classdef SEIRModel
    % # Description
    % This class solves the SEIR epidemiology models. The implementation
    % here is based on: [IDM docs]
    % (https://docs.idmod.org/projects/emod-hiv/en/latest/model-seir.html)
    % Refer to the README for detailed documentation
    % 
    % # License
    % Copyright (C) 2021  Soham Jariwala
    % 
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    %
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License
    % along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
    properties
        %% Parameters
        N;      % Total number of people
        mu;     % Birth/death rate (assumed equal)
        beta;   % Contact rate
        gamma;  % Infection frequency (1/infection period)
        a;      % Latent frequency (1/latent period) 
        xi;     % Lost immunity rate
        R0;     % Reproduction number
        modelType; % Type of model--closed/vital dyanmics SEIR(S)
    end
    
    methods
      %% Generate model populated with default parameters
        function obj = Constructor(obj)
           % Constructor: A class constructor method for assinging params
           %    Assign default values to the parameters that appear in the
           %    SEIR model.
           
           obj.N = 200;
           obj.mu = 0.005;
           obj.gamma = 0.035;
           obj.beta = 0.1;
           obj.a = 0.1;
           obj.xi = 0.001;
           obj.modelType = 'closedSEIR'; % SEIR Model without vital dynamics
        end
      
      %% Evaulate the reproduction number
        function x = get.R0(obj)
            % R0: Calculation of R0 for the given set of parameters
            %   The reproduction number can be determined directly from the
            %   parameters. It can also be used to determine the herd
            %   immunity population.
            
            mu = obj.mu; beta = obj.beta; gamma = obj.gamma; a = obj.a;
            N = obj.N;
            
            switch(obj.modelType)
                case 'closedSEIR'
                    x = beta/gamma;
                    
                case 'vitalSEIR'
                    x = a*beta/(mu + a)/(mu + gamma);
             
                case 'closedSEIRS'
                    x = beta/gamma;
                
                case 'vitalSEIRS'
                    x = a*beta/(mu + a)/(mu + gamma);
            end
        end        
        
      %% Model equations
        function dSEIR = Equations(obj,~,SEIR)
            % Equations: dynamic differential equations for the SEIR model
            %   The differential equations describe the evolution of the
            %   infection with an incubation periodic. Four population
            %   compartments are chosen, namely, Susceptible (S), Exposed
            %   (E), Infected (I), Recovered (R).
            
            % Assigning parameters and variables
            mu = obj.mu; beta = obj.beta; gamma = obj.gamma; a = obj.a;
            N = obj.N; xi = obj.xi;

            S = SEIR(1); E = SEIR(2); I = SEIR(3); R = SEIR(4);
            
            % Model equations
            switch(obj.modelType)
                case 'closedSEIR'
                    dS =  - beta*I*S/N;
                    dE = beta*I*S/N - a*E;
                    dI = a*E - gamma*I;
                    dR = gamma*I;
                
                case 'vitalSEIR'
                    dS = mu*N - mu*S - beta*I*S/N;
                    dE = beta*I*S/N - (mu + a)*E;
                    dI = a*E - (gamma + mu)*I;
                    dR = gamma*I - mu*R;
                
                case 'closedSEIRS'
                    dS = - beta*I*S/N + xi*R;
                    dE = beta*I*S/N - a*E;
                    dI = a*E - gamma*I;
                    dR = gamma*I - xi*R;
                
                case 'vitalSEIRS'
                    dS = mu*N - mu*S - beta*I*S/N + xi*R;
                    dE = beta*I*S/N - (mu + a)*E;
                    dI = a*E - (gamma + mu)*I;
                    dR = gamma*I - mu*R - xi*R;
            end
            
            % Output derivative column vector
            dSEIR = [dS; dE; dI; dR]; 
        end
      
      %% Solution for equilibrium and time dependent (Dynamic) behavior
        function eqSEIR = Equilibrium(obj,currentPopulation)
            % Equilibrium: Solves for the equilibrium condition
            %   The SEIR models reach an endemic equilibrium at R0 > 1 and
            %   a disease free equilibrium at R0 < 1. One can find out by
            %   solving the equation for steady state.
            if ~strcmp(obj.modelType,'vitalSEIR' ) && ~strcmp(obj.modelType,'vitalSEIRS' )
                fprintf("No equilibrium population found: Multiple steady states possible\n\n");
                eqSEIR = [obj.N; 0; 0; 0];
                return
            end
            mu = obj.mu; beta = obj.beta; gamma = obj.gamma; a = obj.a;
            N = obj.N; xi = obj.xi;

            fun = @(seir) obj.Equations(0,seir);
            
            currentPopulation = [0.4; 0.2; 0.2; 0.2]*N;
            
            % Solve for the steady state to get the equilibrium
            options = optimset('TolFun', eps, 'TolX', eps, 'Display','off');
            
                [eqSEIR,~,EXITFLAG,~,~] = fsolve(fun, currentPopulation, options);
        end
      
        function SEIR = Dynamic(obj, time, init)
            % Dynamic: Solves for the time dependent condition
            %   The SEIR model has a time dependent behavior where one can
            %   observe the evolution of Susceptible (S), Exposed
            %   (E), Infected (I), and Recovered (R) populations.
                      
            mu = obj.mu; beta = obj.beta; gamma = obj.gamma; a = obj.a;
            N = obj.N; xi = obj.xi;
            
            if nargin < 3
                fprintf('Initial value not provided: Default initial value used\n\n');
                init = [obj.N-1;1;0;0];
            end 
        
            odeopts = odeset('RelTol',1e-8, 'Stats','off');
            try
            sol = ode15s(@obj.Equations, [time(1), time(end)], init, odeopts);
            SEIR = deval(sol, time)';

            catch
            SEIR = zeros(length(time), 4);
            end
        end
      
      %% Data fitting
        function obj = FitData(obj, optimType, data, init)
           % FitData: Function that takes data and fits a SEIR(S) model
           %    Obtain the model parameters and R0 for a time series data
           %    using a dynamic data fitting. The data needs to be in the
           %    format [time, S, E, I, R], where each entry is a column
           %    vector.
           
           time = data(:,1); S_Exp = data(:,2); E_Exp = data(:,3);
           I_Exp = data(:,4); R_Exp = data(:,5);
           
           function out = FObj(x)
               if strcmp(obj.modelType,'closedSEIR')
                       obj.gamma = x(1);
                       obj.beta  = x(2);
                       obj.a     = x(3);

               elseif strcmp(obj.modelType,'closedSEIRS')
                       obj.gamma = x(1);
                       obj.beta  = x(2);
                       obj.a     = x(3);
                       obj.xi    = x(4);

               elseif strcmp(obj.modelType,'vitalSEIR')
                       obj.gamma = x(1);
                       obj.beta  = x(2);
                       obj.a     = x(3);
                       obj.mu    = x(4);

               elseif strcmp(obj.modelType,'vitalSEIRS')
                       obj.gamma = x(1);
                       obj.beta  = x(2);
                       obj.a     = x(3);
                       obj.mu    = x(4);
                       obj.xi    = x(5);

               end
                
                SEIR = obj.Dynamic(time, init);
                out = sqrt(mean(((S_Exp - SEIR(:,1))./mean(S_Exp)).^2 ...
                               +((E_Exp - SEIR(:,2))./mean(E_Exp)).^2 ...
                               +((I_Exp - SEIR(:,3))./mean(I_Exp)).^2 ...
                               +((R_Exp - SEIR(:,4))./mean(R_Exp)).^2));
           end
           
           % Assigning the intial value to the solution variable 
           if strcmp(obj.modelType,'closedSEIR')
            x0 = zeros(3,1);
                   x0(1) = obj.gamma;
                   x0(2) = obj.beta;
                   x0(3) = obj.a;
                   
           elseif strcmp(obj.modelType,'closedSEIRS')
            x0 = zeros(4,1);
                   x0(1) = obj.gamma;
                   x0(2) = obj.beta;
                   x0(3) = obj.a;
                   x0(4) = obj.xi;
                   
           elseif strcmp(obj.modelType,'vitalSEIR')
            x0 = zeros(4,1);
                   x0(1) = obj.gamma;
                   x0(2) = obj.beta;
                   x0(3) = obj.a;
                   x0(4) = obj.mu;
               
           elseif strcmp(obj.modelType,'vitalSEIRS')
            x0 = zeros(5,1);
                   x0(1) = obj.gamma;
                   x0(2) = obj.beta;
                   x0(3) = obj.a;
                   x0(4) = obj.mu;
                   x0(5) = obj.xi;
                   
           end
           lb = zeros(size(x0));
           ub = ones(size(x0));
           
           if strcmp(optimType, 'partemp')
               fprintf("Using Parallel Tempering Algorithm for Optimization\n\n")
               X = ParTemp(@FObj,x0',...
                   1,...
                   lb','MIN',...
                   ub','MAX');
           elseif strcmp(optimType, 'simanneal') || strcmp(optimType, '')
               fprintf("Using MATLAB:simulannealbnd for Optimization\n\n");
               options = optimoptions('simulannealbnd',...
                   'FunctionTolerance', eps); 
               [X, FEVAL, EXITFLAG] = simulannealbnd(@FObj, ...
                   x0',...
                   lb',ub',...
                   options) %#ok
           end
           
           % Setting the parameters to the best value obtained
           if strcmp(obj.modelType,'closedSEIR')
                   obj.gamma = X(1);
                   obj.beta  = X(2);
                   obj.a     = X(3);

           elseif strcmp(obj.modelType,'closedSEIRS')
                   obj.gamma = X(1);
                   obj.beta  = X(2);
                   obj.a     = X(3);
                   obj.xi    = X(4);

           elseif strcmp(obj.modelType,'vitalSEIR')
                   obj.gamma = X(1);
                   obj.beta  = X(2);
                   obj.a     = X(3);
                   obj.mu    = X(4);

           elseif strcmp(obj.modelType,'vitalSEIRS')
                   obj.gamma = X(1);
                   obj.beta  = X(2);
                   obj.a     = X(3);
                   obj.mu    = X(4);
                   obj.xi    = X(5);

           end
        end
      %% Plot function
        function plotDynamic(obj, time, init)
        % plotDynamic: plots the time dependent condition

            if nargin < 3
                fprintf('Initial value not provided: Default initial value used\n\n');
                init = [obj.N-1;1;0;0];
            end 
            SEIR = Dynamic(obj, time, init);
            
            S = SEIR(:,1); E = SEIR(:,2); I = SEIR(:,3); R = SEIR(:,4);
            
            % Plot dynamic data
            plot(time, S, 'g', ...
                 time, E, 'm',...
                 time, I, 'r',...
                 time, R, 'b',...
                'linewidth',2 ...
                ); hold on;
            
            % Plot Equilibrium values
            if strcmp(obj.modelType,'vitalSEIR' ) || strcmp(obj.modelType,'vitalSEIRS' )
            eqSEIR = Equilibrium(obj, init);
            plot([time(1) time(end)], [eqSEIR(1) eqSEIR(1)],'g-.',...
                 [time(1) time(end)], [eqSEIR(2) eqSEIR(2)],'m-.',...
                 [time(1) time(end)], [eqSEIR(3) eqSEIR(3)],'r-.',...
                 [time(1) time(end)], [eqSEIR(4) eqSEIR(4)],'b-.')
            end
            
            % Figure properties
            set(gca,'FontSize',14,...
                'linewidth',2,...
                'FontName','Times');            
            legend('Susceptible (S)', 'Exposed (E)', 'Infected (I)', ...
                'Recovered (R)', 'location', 'best');
            xlabel('Time (days)');
            ylabel('Number');
            xlim([time(1) time(end)]);
            ylim([0 obj.N]);
        end
    end
end