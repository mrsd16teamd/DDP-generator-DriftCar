function [Fx,Fy] = tire_dyn(Ux, Ux_cmd, mu, mu_slide, Fz, C_x, C_alpha, alpha)

    % longitude wheel slip
%     if (Ux_cmd == Ux)
%         K = 0;
%     elseif Ux == 0
%         Fx = sign(Ux_cmd)*mu*Fz;
%         Fy = 0;
%         return;
%     else
%         K = (Ux_cmd-Ux)/abs(Ux);
%     end
%     

    K = (Ux_cmd-Ux)/(abs(Ux)+1e-3);
    % instead of avoiding -1, now look for positive equivalent
    reverse = 1;
    if K < 0
        reverse = -1;
        K = abs(K);
    end
    
    % alpha > pi/2 cannot be adapted to this formula
    % because of the use of tan(). Use the equivalent angle instead.
    if abs(alpha) > pi/2
        alpha = (pi-abs(alpha))*sign(alpha);
    end
 
    gamma = sqrt(C_x^2*(K/(1+K))^2+C_alpha^2*(tan(alpha)/(1+K))^2);
    
    if gamma <= 3*mu*Fz
        F = 1 - 1/(3*mu*Fz)*(2-mu_slide/mu)*gamma + 1/(9*mu^2*Fz^2)*(1-(2/3)*(mu_slide/mu))*gamma^2;
    else
        % more accurate modeling with peak friction value
        F = mu_slide*Fz/gamma;
    end
    
    Fx = C_x * (K/(1+K)) * F * reverse;
    Fy = -C_alpha * (tan(alpha)/(1+K)) * F;
end