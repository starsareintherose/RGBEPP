import std.stdio;
import std.file;
import std.algorithm;
import std.array;
import std.conv;
import std.getopt;
import std.path;

void main(string[] args)
{
    if (args.length < 4 || args.length > 5 )
    {
        writeln("Usage: "~ args[0] ~" <folder> <condition> <number> [output_file]");
        return;
    }

    // getargs
    string folder = args[1];
    string condition = args[2];
    long threshold = to!long(args[3]);
    bool hasOutputFile = args.length > 4;
    string outputFile = hasOutputFile ? args[4] : "";

    // resotre results
    string result;

    // ie all the file
    foreach (entry; dirEntries(folder, SpanMode.shallow))
    {
        if (entry.isFile)
        {
            // count  '>' 
            size_t count = 0;
            foreach (line; File(entry.name).byLine())
            {
                if (line.length > 0 && line[0] == '>')
                {
                    count++;
                }
            }

            // judge
            bool conditionMet = false;
            if (condition == "=" && count == threshold) conditionMet = true;
            else if (condition == "<" && count < threshold) conditionMet = true;
            else if (condition == ">" && count > threshold) conditionMet = true;
            else if (condition == "<=" && count <= threshold) conditionMet = true;
            else if (condition == ">=" && count >= threshold) conditionMet = true;
	    else if (condition == "><" && count != threshold) conditionMet = true;
            else if (condition == "<>" && count != threshold) conditionMet = true;


            if (conditionMet)
            {
                result ~= baseName(entry.name) ~ " : " ~ to!string(count) ~ "\n";
            }
        }
    }

    // write
    if (hasOutputFile)
    {
	// write to file
        std.file.write(outputFile, result);
    }
    else
    {
        // output to consol
        writeln(result);
    }
}

