function relerr = relative_error(xref, xest)
% Relative error between xref (reference) and xest, up to circular shifts.

    xest = align_to_reference(xest, xref);
    relerr = norm(xref-xest) / norm(xref);

end
