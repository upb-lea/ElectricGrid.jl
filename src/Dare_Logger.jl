using LoggingExtras
using Dates

function DareLogger(;file_name = "", add_date = true, log_file = false, log_level = Logging.Info)

    isdir("log") || mkdir("log")

    filename = "log/"
    filename *= file_name

    if add_date
        filename *= "_" * replace(string(now())[1:end-4], ":" => "-")
    end

    filename *= ".log"

    if log_file
        return TeeLogger(
            MinLevelLogger(FileLogger(filename), log_level),
            MinLevelLogger(ConsoleLogger(), log_level),
        )
    else
        return MinLevelLogger(ConsoleLogger(), log_level)
    end
end