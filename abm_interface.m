function [ interface ] = abm_interface()
    %Interfaces Declaration
    interface = struct(...
        'n1', 1, ...
        'n2', 1, ...
        'perturbanceDownTop', inf, ...
        'perturbanceDownBottom', inf, ...
        'perturbanceUpTop', inf, ...
        'perturbanceUpBottom', inf, ...
        'thickness', 0, ...
        'absorptionCoefficient', 0, ...
        'name', 'unnamed' ...
    );
end
