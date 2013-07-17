
options(
    warn=2, # This will change all warnings into errors, so warnings will also be handled like errors
    error= quote({
      q("no",status=10,FALSE) # standard way for R to end after errors
    })
)

if( runif(1) < 1.5 )	{ stop('Error'); }
cat('OK\n');

