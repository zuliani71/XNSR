function S = add_nested_field(S, fields, value)
    % Aggiunge un campo annidato senza sovrascrivere la struttura esistente
    current_field = fields{1};

    if numel(fields) == 1
        % Caso base: aggiungi il campo finale
        S.(current_field) = value;
    else
        % Se il campo esiste già ed è una struct, continua la ricorsione
        if isfield(S, current_field) && isstruct(S.(current_field))
            S.(current_field) = add_nested_field(S.(current_field), fields(2:end), value);
        else
            % Crea la sottostruttura se non esiste
            S.(current_field) = add_nested_field(struct(), fields(2:end), value);
        end
    end
end