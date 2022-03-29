module FlowSystem
    export Equation, system!

    """
    Custom struct for topography/IC functions
    """
    struct Equation
        f::Function
        params::Tuple
    end

    """
    Builds the system of ODEs (represented by du)
    u = (h_-1, h_0, ..., h_nx, h_nx+1) where h_-1 and h_nx+1 are ghost points
    """
    function system!(du, u, p, t)
        D, alpha, dx, C, k_i, topo, params = p

        function phi(i)
            u[i] + topo(dx*(i-2), params...)
        end
        function interp(i)
            1/2 * (u[i]^3 + u[i+1]^3)
        end
        function f1(i)
            1/(dx^2) * (interp(i-1)*(phi(i-1) - phi(i)) + interp(i)*(phi(i+1) - phi(i)))
        end
        function f2(i)
            (1/dx^4) * (interp(i-1)*(phi(i-2) - 3*phi(i-1) + 3*phi(i) - phi(i+1)) + interp(i)*(-phi(i-1) + 3*phi(i) - 3*phi(i+1) + phi(i+2)))
        end
        function f3(i)
            (1/(2*dx)) * (u[i+1]^3 - u[i-1]^3)
        end
        function f4(i)

            function expo(j)
                x_j = dx*(j-2)
                return exp(2*k_i*x_j)
            end

            (1/(2*dx)) * ((u[i+1]^3 * expo(i+1)) - (u[i-1]^3 * expo(i-1)))
        end

        du[1] = 0
        du[2] = 0
        du[end-1] = 0
        du[end] = 0
        for i in 3:length(u)-2
            #x_i = dx*(i-2)
            #expo = C*exp(2*k_i*x_i)
            du[i] = D*cos(alpha)*f1(i) - f2(i) - sin(alpha)*f3(i) + (2*C*k_i)*f4(i)
        end
        return du
    end
end