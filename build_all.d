#!/usr/bin/env rdmd

import std.process : environment, spawnProcess, wait;
import std.stdio : stderr;
import std.exception : enforce;
import std.array : array;

int main(string[] args)
{
    const compiler = environment.get("DUB_COMPILER", "ldc2");
    const buildType = environment.get("DUB_BUILD_TYPE", "release-ldc");

    immutable string[] configs = [
        "rgbepp",
        "rgbepp-refmix",
        "sortdiamond",
        "splitfasta",
        "countTaxa",
        "delstop",
        "deltaxa",
        "concataln"
    ];
 
    const extraArgs = args[1 .. $];

    foreach (config; configs)
    {
        string[] command = [
            "dub",
            "build",
            "--compiler=" ~ compiler,
            "--build=" ~ buildType,
            "--config=" ~ config
        ];

        command ~= extraArgs;

        stderr.writefln(
            "Building configuration '%s' with compiler '%s' and build type '%s'...",
            config,
            compiler,
            buildType
        );

        auto process = spawnProcess(command);
        const status = wait(process);

        if (status != 0)
        {
            stderr.writefln(
                "Build failed for configuration '%s' with exit code %s.",
                config,
                status
            );
            return status;
        }
    }

    return 0;
}
