using LoggingExtras
using Dates

function DareLogger(;file_name = "", add_timestamp = true, log_file = false, log_level = Logging.Info, log_level_file = Logging.Debug)

    if log_file
        isdir("log") || mkdir("log")
    end

    filename = "log/"
    filename *= file_name

    if add_timestamp
        filename *= "_" * replace(string(now())[1:end-4], ":" => "-")
    end

    filename *= ".log"

    if isnothing(log_level_file)
        log_level_file = log_level
    end

    if log_file
        return TeeLogger(
            MinLevelLogger(FileLogger(filename), log_level_file),
            MinLevelLogger(ConsoleLogger(), log_level),
        )
    else
        return MinLevelLogger(ConsoleLogger(), log_level)
    end
end