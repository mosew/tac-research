function fk = fk_from_ak_efs(ak,efs)
    fk = @(s) sum(ak.*cellfun(@(c) feval(c,s),efs));
end