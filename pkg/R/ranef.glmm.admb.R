"ranef.glmm.admb" <-
function (object, sd = FALSE)
{
    if (!sd)
        out <- object$U
    else out <- object$sd_U
    out
}
