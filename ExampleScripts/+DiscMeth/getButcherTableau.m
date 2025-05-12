%==========================================================================
%> @file getButcherTableau.m
%> @brief Returns Butcher tableau for selected implicit Runge-Kutta method
%>
%> @details
%> This function provides the Butcher tableau (A, b, c) for common
%> Runge-Kutta methods. The tableau defines the coefficients used
%> in the time integration scheme.
%>
%> @param method    String: Name of the Runge-Kutta method
%>
%> @retval A        Matrix: Coefficient matrix A of the Butcher tableau
%> @retval b        Vector: Weight vector b of the Butcher tableau
%> @retval c        Vector: Node vector c of the Butcher tableau
%>
%> @author  Alexander Weiss
%> @date    May, 2025
%==========================================================================
function [A, b, c, name] = getButcherTableau(discMeth)
    switch discMeth
        case 'BE'
            % Backward Euler (s=1, p=1)
            A = 1;
            b = 1;
            c = 1;
            name = 'Backward Euler';
    
        case 'ME'
            % Midpoint Method (s=1, p=2)
            A = 1/2;
            b = 1;
            c = 1/2;
            name = 'Midpoint Euler';

        case 'LIIIc-2'
            % Lobatto IIIc (s=2, p=2) 
            A = [1/2, -1/2;
                1/2, 1/2];
            b = [1/2; 1/2];
            c = [0; 1];
            name = 'Lobatto IIIc with two stages';


        case 'RIIa-2'
            % Radau IIa (s=2, p=3) 
            A = [5/12, -1/12;
                3/4, 1/4];
            b = [3/4; 1/4];
            c = [1/3; 1];
            name = 'Radau IIa with two stages';

    
        case 'LIIIc-3'
            % Lobatto IIIc (s=3, p=4)
            A = [1/6, -1/3, 1/6;
                 1/6, 5/12, -1/12;
                 1/6, 2/3, 1/6];
            b = [1/6; 2/3; 1/6];
            c = [0; 1/2; 1];
            name = 'Lobatto IIIc with three stages';

        case 'RIIa-3'
            % Radau IIa (s=2, p=5)
            sqrt6 = sqrt(6);

            A = [ (88 - 7*sqrt6)/360,    (296 - 169*sqrt6)/1800, (-2 + 3*sqrt6)/225;
                (296 + 169*sqrt6)/1800,  (88 + 7*sqrt6)/360,     (-2 - 3*sqrt6)/225;
                (16 - sqrt6)/36,         (16 + sqrt6)/36,         1/9 ];

            b = [(16 - sqrt6)/36;        (16 + sqrt6)/36;       1/9] ;  

            c = [(4 - sqrt6)/10;
                (4 + sqrt6)/10;
                1];
            name = 'Radau IIa with three stages';


        otherwise
            error('Unknown method "%s".', discMeth);
    end
end